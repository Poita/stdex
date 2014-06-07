module stdex.statistics;

import std.algorithm;
import std.conv;
import std.math;
import std.mathspecial;
import std.random;
import std.range;
import stdex.algorithm;
import stdex.math;
import stdex.range;
import stdex.util;

version(unittest)
{
    import std.stdio;
}

/**
Computes the arithmetic mean of the elements of $(D samples).
This is equal to the sum of the elements divided by the number of elements.
$(D samples) must be a non-empty, finite, input range.
 */
auto mean(Range)(Range samples)
{
    // TODO: Use Kahan.

    static assert(isInputRange!Range, "samples must be an input range.");
    static assert(!isInfinite!Range, "samples must not be an infinite range.");

    assert(!samples.empty, "Cannot compute mean of an empty range.");
    auto sum = samples.front;
    size_t n = 1;
    samples.popFront();
    while (!samples.empty)
    {
        sum += samples.front;
        ++n;
        samples.popFront();
    }
    return sum / n;
}

///
unittest
{
    assert(mean([1, 3, 2]) == 2);
}

/**
Computes the variance of the $(D samples) from the $(D mean).
By default, $(D variance) computes the biased sample variance.
To compute the unbiased sample variance, set $(D sampleCorrection) to $(D 1).
$(D samples) must be a non-empty, finite, input range.
 */
real variance(Range, T)(Range samples, T mean, size_t sampleCorrection = 0)
{
    // TODO: Use better numerical algo.
    // See: http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    static assert(isInputRange!Range, "samples must be an input range.");
    static assert(!isInfinite!Range, "samples must not be an infinite range.");

    assert(!samples.empty, "Cannot compute standard deviation of an empty range.");
    assert(sampleCorrection <= 1, "sampleCorrection must be either 0 or 1");

    real sum = sqr(samples.front - mean);
    size_t n = 1 - sampleCorrection;
    samples.popFront();
    while (!samples.empty)
    {
        sum += sqr(samples.front - mean);
        ++n;
        samples.popFront();
    }
    return sum / n;
}

/**
Computes the standard deviation of the $(D samples) from the $(D mean).
This is equal to the square root of the variance.
By default, $(D standardDeviation) computes the biased sample standard deviation.
To compute the unbiased sample standard deviation, set $(D sampleCorrection) to $(D 1).
$(D samples) must be a non-empty, finite, input range.
 */
real standardDeviation(Range, T)(Range samples, T mean, size_t sampleCorrection = 0)
{
    return sqrt(variance(samples, mean, sampleCorrection));
}

/**
Computes the median of the elements of $(D samples).
When there is more than one median, any may be returned.
$(D samples) must be a non-empty, finite, input range.
 */
T median(Range)(Range samples)
{
    // TODO: should use a predicate for sorting.
    static assert(isInputRange!Range, "samples must be an input range.");
    static assert(!isInfinite!Range, "samples must not be an infinite range.");
    assert(!samples.empty, "Cannot compute the median of an empty sample.");

    // TODO: use the O(n) selection algorithm
    return samples.array.sort()[$ / 2];
}

/**
Computes the mode of the elements of $(D samples).
When there is more than one mode, any may be returned.
$(D samples) must be a non-empty, finite, input range.
 */
ElementType!Range mode(Range)(Range samples)
{
    // TODO: should use a predicate for sorting.
    static assert(isInputRange!Range, "samples must be an input range.");
    static assert(!isInfinite!Range, "samples must not be an infinite range.");
    assert(!samples.empty, "Cannot compute the mode of an empty sample.");

    return samples.array.sort().group.argMax!("a[1]")[0];
}

///
unittest
{
    assert([0, 2, 1, 1].mode == 1);
    assert([0].mode == 0);
    assert([1, 1, 0, 0, 0, 2].mode == 0);
}

/**
Computes the unbiased sample standard deviation.
$(D samples) must be a non-empty, finite, input range.
 */
real sampleStandardDeviation(Range, T)(Range samples, T mean)
{
    return standardDeviation(samples, mean, 1);
}

/**
Computes the standard error of a sampling distribution of size $(D sampleSize)
and standard deviation $(D stdDeviation).
 */
real standardError(real stdDeviation, size_t sampleSize)
{
    assert(sampleSize != 0, "Cannot compute the standard error of an empty sampling distribution.");
    return stdDeviation / sqrt(sampleSize.to!real);
}

/**
Computes the p-value for an observation. TODO
 */
real pValue(real populationMean, real sampleMean, real sampleStdDeviation, size_t sampleSize)
{
    real error = standardError(sampleStdDeviation, sampleSize);
    auto difference = abs(populationMean - sampleMean);
    if (sampleSize > 30)
    {
        // Large sample size, use Normal Distribution.
        auto samplingDist = NormalDistribution(0.0, 1.0);
        return 2 * samplingDist.cdf(-difference / error);
    }
    else
    {
        // Small sample size, use Student's T-Distribution
        const degreesOfFreedom = sampleSize - 1;
        auto samplingDist = StudentTDistribution(degreesOfFreedom);
        return 2 * samplingDist.cdf(-difference / error);
    }
}

auto chiSquaredTest(ObservedRange, ExpectedRange)
                   (ObservedRange observed, ExpectedRange expected, size_t degreesOfFreedom = size_t.max)
{
    real x2 = 0.0;
    size_t calculatedDoF = 0;
    foreach (e; zip(observed, expected))
    {
        x2 += sqr(e[0] - e[1]) / e[1];
        ++calculatedDoF;
    }
    size_t dof = degreesOfFreedom == size_t.max ? calculatedDoF : degreesOfFreedom;
    real p = ChiSquaredDistribution(dof).sf(x2);
    return namedTuple!("stat", "p")(x2, p);
}

auto oneWayFRatio(RangeOfSamples)(RangeOfSamples sampleSets)
{
    // TODO: could be done faster, but this is easier to read and fast enough.
    real[] groupMeans = sampleSets.save.map!mean.array();
    size_t[] groupSizes = sampleSets.save.map!walkLength.array();
    real overallMean = groupMeans.mean();
    real sB = zip(groupMeans, groupSizes).map!(x => x[1] * sqr(x[0] - overallMean)).sum();
    real sW = zip(sampleSets, groupMeans).map!(x => x[0].map!(s => sqr(s - x[1])).sum()).sum();
    size_t fB = sampleSets.walkLength - 1;
    size_t fW = groupSizes.map!(x => x - 1).sum();
    real fRatio = (sB * fW) / (sW * fB);
    real p = FisherFDistribution(fB, fW).sf(fRatio);
    return namedTuple!("stat", "p")(fRatio, p);
}

///
unittest
{
    real[][] samples = [[6, 8, 4, 5, 3, 4],
                        [8, 12, 9, 11, 6, 8],
                        [13, 9, 11, 8, 7, 12]];
    auto r = oneWayFRatio(samples);
    assert(r.stat.approxEqual(9.26471));
    assert(r.p.approxEqual(0.00239878));
}

/**
Models a continuous uniform probability distribution, with constant density
over the [$(D min), $(D max)) range.
 */
struct ContinuousUniformDistribution
{
public:
    this(real min = 0, real max = 1)
    {
        assert(min <= max, "min must be less than or equal to max.");
        m_min = min;
        m_max = max;
    }

    real pdf(real x)
    {
        return x < m_min ? 0.0 :
               x >= m_max ? 0.0 :
               1.0 / (m_max - m_min);
    }

    real cdf(real x)
    {
        return inverseLerpClamp(x, m_min, m_max);
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        return uniform(m_min, m_max, urng);
    }

    @property real min()
    {
        return m_min;
    }

    @property real max()
    {
        return m_max;
    }

    @property real mean()
    {
        return 0.5 * (m_min + m_max);
    }

    @property real mode()
    {
        return mean;
    }

    @property real median()
    {
        return mean;
    }

    @property real stdDeviation()
    {
        enum real c_sqrt12 = sqrt(12.0);
        return (m_max - m_min) * c_sqrt12;
    }

    @property real variance()
    {
        return sqr(m_max - m_min) / 12.0;
    }

    @property real skewness()
    {
        return 0;
    }

    @property real kurtosis()
    {
        enum real c_kurtosis = -6.0 / 5.0;
        return c_kurtosis;
    }

    @property real entropy()
    {
        return log(m_max - m_min);
    }

private:
    real m_min;
    real m_max;
}

struct NormalDistribution
{
public:
    this(real mean, real stdDeviation)
    {
        m_mean = mean;
        m_stdDeviation = stdDeviation;
    }

    real pdf(real x)
    {
        enum real factor = 1 / sqrt(2 * PI);
        x -= m_mean;
        return factor / (m_stdDeviation * exp(x * x / (2 * variance)));
    }

    real cdf(real x)
    {
        return std.mathspecial.normalDistribution((x - m_mean) / m_stdDeviation);
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        real ret;
        if (m_nextSample == m_nextSample)
        {
            ret = m_nextSample;
            m_nextSample = real.nan;
        }
        else
        {
            real x, y, r;
            do
            {
                x = uniform(-1.0, 1.0, urng);
                y = uniform(-1.0, 1.0, urng);
                r = x * x + y * y;
            }
            while (r > 1.0 || r == 0.0);

            const m = sqrt(-2.0 * log(r) / r);
            m_nextSample = x * m;
            ret = y * m;
        }
        ret *= m_stdDeviation;
        ret += m_mean;
        return ret;
    }

    @property real mean()
    {
        return m_mean;
    }

    @property real mode()
    {
        return m_mean;
    }

    @property real median()
    {
        return m_mean;
    }

    @property real stdDeviation()
    {
        return m_stdDeviation;
    }

    @property real variance()
    {
        return m_stdDeviation * m_stdDeviation;
    }

    @property real skewness()
    {
        return 0;
    }

    @property real kurtosis()
    {
        return 0;
    }

    @property real entropy()
    {
        return 0.5 * log(2 * PI * E * variance);
    }

private:
    real m_mean;
    real m_stdDeviation;
    real m_nextSample = real.nan;
}

struct StudentTDistribution
{
public:
    this(real degreesOfFreedom = 1)
    {
        m_nu = degreesOfFreedom;
        m_normalD = NormalDistribution(0.0f, 1.0);
        m_gammaD = GammaDistribution(m_nu / 2.0, 2.0);
    }

    real pdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        real v = (m_nu + 1.0) / 2.0;
        real f = gamma(v);
        f /= sqrt(m_nu * PI) * gamma(m_nu / 2.0);
        f *= pow(1.0 + sqr(x) / m_nu, -v);
        return f;
    }

    real cdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        // TODO: write hypergeometric
        real halfNu = m_nu / 2.0;
        real c = 1.0 - 0.5 * betaIncomplete(halfNu, 0.5, m_nu / (sqr(x) + m_nu));
        return x > 0 ? c : 1 - c; 
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        return m_normalD.sample(urng) * sqrt(m_nu / m_gammaD.sample(urng));
    }

    @property real nu()
    {
        return m_nu;
    }

    @property real mean()
    {
        return m_nu > 1 ? 0 :
               real.nan;
    }
    
    @property real mode()
    {
        return 0;
    }
    
    @property real median()
    {
        return 0;
    }
    
    @property real stdDeviation()
    {
        return sqrt(variance);
    }

    @property real variance()
    {
        return m_nu > 2 ? m_nu / (m_nu - 2) :
               m_nu > 1 ? real.infinity :
               real.nan;
    }

    @property real skewness()
    {
        return m_nu > 3 ? 0 :
               real.nan;
    }
    
    @property real kurtosis()
    {
        return m_nu > 4 ? 6 / (m_nu - 4) :
               m_nu > 2 ? real.infinity :
               real.nan;
    }

    @property real entropy()
    {
        real e = (m_nu + 1) / 2;
        e *= digamma(e) - digamma(m_nu / 2);
        e += log(sqrt(m_nu) * beta(m_nu / 2, 0.5));
        return e;
    }

private:
    real m_nu;
    NormalDistribution m_normalD;
    GammaDistribution m_gammaD;
}

struct ChiSquaredDistribution
{
public:
    this(real degreesOfFreedom = 1)
    {
        m_k = degreesOfFreedom;
        m_gammaD = GammaDistribution(degreesOfFreedom / 2.0, 1.0);
    }

    real pdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        real halfK = k / 2.0;
        return pow(x, halfK - 1.0) * exp(-x / 2.0) / pow(2.0, halfK) * gamma(halfK);
    }

    real cdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        real halfK = k / 2.0;
        return gammaIncomplete(halfK, x / 2.0) / gamma(halfK);
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        return 2.0 * m_gammaD.sample(urng);
    }

    @property real k()
    {
        return m_k;
    }

    @property real mean()
    {
        return m_k;
    }

    @property real mode()
    {
        return max(k - 2.0, 0.0);
    }

    @property real median()
    {
        // TODO: this is approximate, is there a better formula?
        return k * pow(1.0 - 2.0 / (9.0 * m_k), 3.0);
    }
    
    @property real stdDeviation()
    {
        return sqrt(variance);
    }

    @property real variance()
    {
        return 2.0 * m_k;
    }

    @property real skewness()
    {
        return sqrt(8.0 / m_k);
    }

    @property real kurtosis()
    {
        return sqrt(12.0 / m_k);
    }

    @property real entropy()
    {
        real halfK = m_k / 2.0;
        return halfK + log(2.0 * gamma(halfK)) + (1.0 - halfK) * digamma(halfK);
    }

private:
    real m_k;
    GammaDistribution m_gammaD;
}

struct FisherFDistribution
{
public:
    this(real d1 = 1, real d2 = 1)
    {
        m_d1 = d1;
        m_d2 = d2;
        m_gammaD1 = GammaDistribution(d1 / 2.0, 1.0);
        m_gammaD2 = GammaDistribution(d2 / 2.0, 1.0);
    }

    real pdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        real p = pow(m_d1 * x, m_d1) * pow(m_d2, m_d2);;
        p /= pow(m_d1 * x + m_d2, m_d1 + m_d2);
        p = sqrt(p);
        p /= x * beta(m_d1 / 2.0, m_d2 / 2.0);
        return p;
    }

    real cdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        return betaIncomplete(m_d1 / 2.0, m_d2 / 2.0, m_d1 * x / (m_d1 * x + m_d2));
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        return (m_gammaD1.sample(urng) * m_d1) / (m_gammaD2.sample(urng) * m_d2);
    }

    @property real m()
    {
        return m_d1;
    }

    @property real n()
    {
        return m_d2;
    }

    @property real mean()
    {
        // TODO: what is this when d2 <= 2.0 ???
        return m_d2 > 2.0 ? m_d2 / (m_d2 - 2.0) : real.nan;
    }

    @property real mode()
    {
        // TODO: what is this when d1 <= 2.0 ???
        return m_d1 > 2.0 ? (m_d1 - 2.0) / m_d1 * (m_d2 / (m_d2 + 2.0)) : real.nan;
    }

    @property real median()
    {
        // TODO:
        return real.nan;
    }
    
    @property real stdDeviation()
    {
        return sqrt(variance);
    }

    @property real variance()
    {
        // TODO: what is this when d2 <= 4.0 ???
        return m_d2 > 4.0 ? (2.0 * sqr(m_d2) * (m_d1 + m_d2 - 2.0)) / (m_d1 * sqr(m_d2 - 2.0) * (m_d2 - 4.0)) : real.nan;
    }

    @property real skewness()
    {
        // TODO: what is this when d2 <= 6.0
        return m_d2 > 6.0 ? ((2.0 * m_d1 + m_d2 - 2.0) * sqrt(8.0 * (m_d2 - 4.0))) /
                            ((m_d2 - 6.0) * sqrt(m_d1 * (m_d1 + m_d2 - 2.0))) : real.nan;
    }

    @property real kurtosis()
    {
        // TODO: what is this when d2 <= 8.0
        real k = m_d1 * (5.0 * m_d2 - 22.0) * (m_d1 + m_d2 - 2.0);
        k += (m_d2 - 4.0) * sqr(m_d2 - 2.0);
        k /= m_d1 * (m_d2 - 6.0) * (m_d2 - 8.0) * (m_d1 + m_d2 - 2.0);
        return 12.0 * k;
    }

    @property real entropy()
    {
        // TODO:
        return real.nan;
    }

private:
    real m_d1;
    real m_d2;
    GammaDistribution m_gammaD1;
    GammaDistribution m_gammaD2;
}

struct GammaDistribution
{
public:
    this(real alpha, real beta)
    {
        assert(alpha > 0, "alpha parameter must be positive.");
        assert(beta > 0, "beta parameter must be positive.");

        m_alpha = alpha;
        m_beta = beta;
        m_malpha = alpha < 1.0 ? alpha + 1.0 : alpha;
        m_alpha2 = 1.0 / sqrt(9.0 * (m_malpha - 1.0 / 3.0));
        m_normalD = NormalDistribution(0.0, 1.0);
    }

    real pdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        assert(x >= 0, "Gamma distribution is not defined for negative numbers.");
        return pow(m_beta, m_alpha) / gamma(m_alpha) * 
               pow(x, m_alpha - 1.0) * exp(-m_beta * x);
    }

    real cdf(real x)
    {
        // TODO: check accuracy of this
        // TODO: cache constants
        assert(x >= 0, "Gamma distribution is not defined for negative numbers.");
        return 1.0 / gamma(m_alpha) * gammaIncomplete(m_alpha, m_beta * x);
    }

    real sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        // Marsaglia, G. and Tsang, W. W.
        // A Simple Method for Generating Gamma Variables
        // ACM Transactions on Mathematical Software, 26, 3, 363-372, 2000.

        real a1 = m_malpha - 1.0 / 3.0;
        real u, v, n, n2;
        do
        {
            do
            {
                n = m_normalD.sample(urng);
                v = 1.0 + m_alpha2 * n;
            }
            while (v <= 0.0);

            v = v * sqr(v);
            u = uniform(0.0, 1.0, urng);
            n2 = sqr(n);
        }
        while(u > 1.0 - 0.331 * sqr(n2) && log(u) > 0.5 * n2 + a1 * (1.0 - v + log(v))); 

        if (m_alpha == m_malpha)
            return a1 * v * m_beta;

        u = uniform!"()"(0.0, 1.0, urng);
        return pow(u, 1.0 / m_alpha) * a1 * v * m_beta;
    }

    @property real alpha()
    {
        return m_alpha;
    }

    @property real beta()
    {
        return m_beta;
    }

    @property real mean()
    {
        return m_alpha / m_beta;
    }
    
    @property real mode()
    {
        // TODO: find a way of calculating this
        return real.nan;
    }
    
    @property real median()
    {
        // TOOD: should do something when alpha <= 1
        return m_alpha > 1 ? (m_alpha - 1) / m_beta : real.nan;
    }
    
    @property real stdDeviation()
    {
        return sqrt(m_alpha) / m_beta;
    }

    @property real variance()
    {
        return m_alpha / sqr(m_beta);
    }

    @property real skewness()
    {
        return 2.0 / sqrt(m_alpha);
    }
    
    @property real kurtosis()
    {
        return 6.0 / m_alpha;
    }

    @property real entropy()
    {
        real e = m_alpha;
        e -= log(m_beta);
        e += logGamma(m_alpha);
        e += (1.0 - m_alpha) * digamma(m_alpha);
        return e;
    }

private:
    real m_alpha;
    real m_beta;
    real m_malpha;
    real m_alpha2;
    NormalDistribution m_normalD;
}


struct BernoulliDistribution
{
public:
    this(real p)
    {
        assert(p >= 0.0 && p <= 1.0, "Bernoulli success probability must be between 0 and 1.");
        m_q = 1 - p;
    }

    real pdf(bool x)
    {
        return x ? 1 - m_q : m_q;
    }

    real cdf(bool x)
    {
        return x ? 1 : m_q;
    }

    bool sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
        return uniform(0.0, 1.0, urng) > m_q;
    }

    @property real p()
    {
        return 1 - m_q;
    }

    @property real q()
    {
        return m_q;
    }

    @property real mean()
    {
        return 1 - m_q;
    }

    @property bool mode()
    {
        // TODO: if m_p == 0.5 there are two modes
        return m_q < 0.5;
    }

    @property real median()
    {
        return m_q > 0.5 ? 0.0 : m_q < 0.5 ? 1.0 : 0.5;
    }

    @property real stdDeviation()
    {
        return sqrt(variance);
    }

    @property real variance()
    {
        return m_q * (1 - m_q);
    }

    @property real skewness()
    {
        return (m_q + m_q - 1.0) / sqrt((1.0 - m_q) * m_q);
    }

    @property real kurtosis()
    {
        return (1.0 - 6.0 * (1.0 - m_q) * m_q) / ((1.0 - m_q) * m_q);
    }

    @property real entropy()
    {
        return -m_q * log(m_q) - (1.0 - m_q) * log(1.0 - m_q);
    }

private:
    real m_q;
}

/+struct BinomialDistribution
{
public:
    this(size_t n, real p)
    {
        assert(p >= 0.0 && p <= 1.0, "Binomial success probability must be between 0 and 1.");
        m_n = n;
        m_p = p;
    }

    real pdf(size_t x)
    {
    }

    real cdf(bool x)
    {
    }

    size_t sample(UniformRandomNumberGenerator)(ref UniformRandomNumberGenerator urng)
    {
    }

    @property real p() { return m_p; }
    @property real q() { return 1 - m_p; }
    @property real mean() { return m_n * m_p; }
    @property size_t mode() { return ; }
    @property real median() { return ; }
    @property real stdDeviation() { return ; }
    @property real variance() { return ; }
    @property real skewness() { return ; }
    @property real kurtosis() { return ; }
    @property real entropy() { return ; }
private:
    size_t m_n;
    real m_p;
}+/

auto sample(Distribution)(Distribution d)
{
    return d.sample(rndGen);
}

auto samples(ProbabilityDistribution, UniformRandomNumberGenerator)
            (ProbabilityDistribution distribution, ref UniformRandomNumberGenerator urng)
{
    return generate!(() => distribution.sample(urng));
}

auto samples(ProbabilityDistribution)(ProbabilityDistribution distribution)
{
    return samples(distribution, rndGen);
}

@property real sf(Distribution)(Distribution dist, real x)
{
    return 1 - dist.cdf(x);
}

/+unittest
{
    import std.stdio;
    writeln("testing stdex.statistics.StudentTDistribution");

    //auto xs = StudentTDistribution(1).samples.take(10000);
    auto xs = ChiSquaredDistribution(3).samples.take(100000);
    int[int] h;
    foreach (x; xs)
        h[cast(int)(floor(x*10))]++;
    auto ha = h.keyValueArray.sort!("a[0] < b[0]");
    foreach (t; ha.filter!("a[1] > 10"))
        writefln("% 3.1f: %s", t[0]/10.0, repeat("*").take(t[1]/50).join);
}+/

// TODO Binomial Distribution
// TODO Categorical Distribution
// TODO Geometric Distribution
// TODO DiscreteUniformDistribution

// Random variates
// Percent point function (inverse of cdf)
// Inverse survivial function (inverse of sf)
// Moment - non-central moments