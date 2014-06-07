module stdex.learning.svm;

import std.algorithm;
import std.math;
import std.random;
import std.range;
import stdex.algorithm;
import stdex.math;
import stdex.util;


class LinearSVM(T)
{
public:
    this(size_t dimensions)
    {
        w.length = dimensions;
    }

    this(T[] w, T b)
    {
        this.w = w.dup;
        this.b = b;
    }

    void updateW(RangeA, RangeX, RangeY)(RangeA a, RangeX x, RangeY y)
    {
        w[] = 0;
        while (!a.empty && !y.empty)
        {
            auto ay = a.front * y.front;
            foreach (k, v; x.front)
                w[k] += ay * v;
            a.popFront();
            x.popFront();
            y.popFront();
        }
    }

    T test(Vector)(in Vector t)
    {
        return dot(w, t) - b;
    }

    @property T[] weights()
    {
        return w;
    }

    @property const T offset()
    {
        return b;
    }

private:
    T[] w;
    T b = 0;
}

class SMO(alias K = dot, Vector, T)
{
public:
    // parameters
    enum T eps = 1e-3; // stopping tolerance
    enum T tau = 1e-12;

    this(Vector[] x, T[] y, T c, size_t dimensions)
    {
        this.svm = new LinearSVM!(T)(dimensions);
        this.x = x;
        this.y = y;
        this.c = c;
        this.a = new T[y.length];
        this.g = new T[y.length];
        this.active = new size_t[y.length];
        this.nactive = 0;
        a[] = 0;
        fillKernelCache();
    }

    typeof(svm) train()
    {
        size_t n = x.length;
        assertf(y.length == n, "Size of x and y must be equal (x.length = %s, y.length = %s)", n, y.length);

        // main routine
        a[] = 0;
        g[] = -1;
        copy(iota(n), active);
        nactive = n;

        int pass = 0;
        import std.stdio;
        while (pass < 1_000_000)
        {
            ++pass;
            if (pass % min(n, 100) == 0)
                shrink();

            size_t i, j;
            selectB(i, j);
            if (j == -1)
            {
                /+g[] = -1;
                foreach (gi; 0..n)
                    foreach (gj; 0..n)
                        g[gi] += a[gj] * kernel(gi, gj);

                ij = selectB();
                i = ij[0];
                j = ij[1];
                if (j == -1)+/
                {
                    break;
                }
            }

            T oldAi = a[i];
            T oldAj = a[j];

            T p = max(kernel(i, i) + kernel(j, j) - 2 * kernel(i, j), tau);
            T q = -y[i] * g[i] + y[j] * g[j];
            T s = y[i] * a[i] + y[j] * a[j];

            a[i] = clamp!T(a[i] + y[i] * q / p, 0, c);
            a[j] = clamp!T(y[j] * (s - y[i] * a[i]), 0, c);
            a[i] = y[i] * (s - y[j] * a[j]);

            T deltaAi = a[i] - oldAi;
            T deltaAj = a[j] - oldAj;
            
            foreach (t; activeSet)
                g[t] += y[i] * y[t] * kernel(t, i) * deltaAi + y[j] * y[t] * kernel(t, j) * deltaAj;
        }

        svm.updateW(a, x, y);
        svm.b = 0;
        //svm.rho = 0;
        size_t nsv = 0;
        foreach (i; 0..n)
            if (a[i] > 0 && a[i] < c)
            {
                ++nsv;
                svm.b += dot(svm.w, x[i]) - y[i];
                //svm.rho += y[i] * g[i];
            }
        if (nsv != 0)
        {
            svm.b /= nsv;
            //svm.rho /= nsv;
        }
        //svm.a = a;
        return svm;
    }

    void selectB(ref size_t outi, ref size_t outj)
    {

        size_t i = activeSet
                   .filter!(k => y[k] == 1 && a[k] < c || y[k] == -1 && a[k] > 0)
                   .argMax!(k => -y[k] * g[k]);
        T gmax = -y[i] * g[i];

        T gmin = 1e20;
        size_t j = -1;
        if (i != -1)
        {
            T obj_min = 1e20;
            T kii = kernel(i, i);
            foreach (t; activeSet)
            {
                if (y[t] == +1 && a[t] > 0 || y[t] == -1 && a[t] < c)
                {
                    T q = gmax + y[t] * g[t];
                    gmin = min(gmin, -y[t] * g[t]);

                    if (q > 0)
                    {
                        T p = max(kii + kernel(t, t) - 2 * kernel(i, t), tau);
                        T s = -(q * q) / p;
                        if (s <= obj_min)
                        {
                            j = t;
                            obj_min = s;
                        }
                    }
                }
            }
        }

        //import std.stdio;
        //writeln(gmax - gmin);
        if (gmax - gmin < eps)
        {
            outi = -1;
            outj = -1;
        }
        else
        {
            outi = i;
            outj = j;
        }
    }

    bool shrinkOne(size_t i, T gmax1, T gmax2)
    {
        if (a[i] == c)
        {
            if(y[i] > 0)
                return -g[i] > gmax1;
            else
                return -g[i] > gmax2;
        }
        else if (a[i] == 0)
        {
            if (y[i] > 0)
                return g[i] > gmax2;
            else        
                return g[i] > gmax1;
        }
        return false;
    }

    void shrink()
    {
        T gmax1 = -1e20;
        T gmax2 = -1e20;

        foreach (i; activeSet)
        {
            if (y[i] == 1)
            {
                if (a[i] != c && -g[i] >= gmax1)
                    gmax1 = -g[i];
                if (a[i] != 0 && g[i] >= gmax2)
                    gmax2 = g[i];
            }
            else
            {
                if (a[i] != c && -g[i] >= gmax2)
                    gmax2 = -g[i];
                if (a[i] != 0 && g[i] >= gmax1)
                    gmax1 = g[i];
            }
        }

        // TODO: unshrink reconstruct gradient?

        size_t k = 0;
        foreach (i; activeSet)
        {
            if (shrinkOne(i, gmax1, gmax2))
            {
                nactive--;
            }
            else
            {
                active[k++] = i;
            }
        }
    }

    auto activeSet() { return active[0..nactive]; }

    auto unboundIndices()
    {
        return iota(y.length).filter!(i => !isBound(i));
    }

    auto supportVectorIndices()
    {
        return iota(y.length).filter!(i => a[i] != 0);
    }

    bool isBound(size_t i)
    {
        return a[i] == 0 || a[i] == c;
    }

    T kernel(size_t i, size_t j)
    {
        size_t idx = i * x.length + j;
        if (kcache[idx] == kcacheSentinel)
            kcache[idx] = K(x[i], x[j]);
        return kcache[idx];
        //return K(x[i], x[j]);
    }

    void fillKernelCache()
    {
        kcache = new T[x.length * x.length];
        kcache[] = kcacheSentinel;
    }

private:
    enum T kcacheSentinel = cast(T)-999;

    LinearSVM!(T) svm;
    T[] a;          // Lagrange multipliers
    Vector[] x;     // Training data vectors
    T[] y;          // Training data labels (-1, +1)
    T c;            // Box conditions (0 <= a[i] <= c)
    T[] g;
    T[] kcache;     // Cache of K(i, j)
    size_t[] active;
    size_t nactive;
}

auto smoTrain(alias K = dot, Vector, T)(Vector[] x, T[] y, T c, size_t dimensions)
{
    return (new SMO!(K, Vector, T)(x, y, c, dimensions)).train();
}

import stdex.profile;

version(none) void main()
{
    import std.stdio;
    writeln("testing stdex.learning.svm");

    //float[][] x = [[1.0f, 0.0f], [0.0f, 1.0f], [2.0f, -1.0f]];
    //float[] y = [-1.0f, 1.0f, -1.0f];
    //float[][] x = [[1.0f, 0.0f], [0.0f, 2.0f], [3.0f, 1.0f], [10.0f, -10.0f], [-10.0f, 20.0f]];
    //float[] y = [-1.0f, 1.0f, -1.0f, -1.0f, 1.0f];

    rndGen.seed(0);

    size_t n = 400;
    double[size_t][] x;
    x.length = n;
    double[] y;
    y.length = n;

    foreach (k; 0..20)
    {
        foreach (i; 0..n)
        {
            y[i] = (i & 1) ? 1.0f : -1.0f;
            if (i & 1)
                x[i] = [0:uniform(-10.0f, 10.0f), 5:uniform(-3.0f, 10.0f)];
            else
                x[i] = [0:uniform(-10.0f, 10.0f), 5:uniform(-10.0f, 3.0f)];
        }
        
        auto svm = smoTrain(x, y, 2.0, 6);
        writefln("w = %s b = %s", svm.w, svm.b);
    }
    version(profile) stdex.profile.flush(File("profile", "w"));
}