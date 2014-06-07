module stdex.learning.decisiontree;

import std.algorithm;
import std.array;
import std.math;
import std.range;
import std.typecons;

real entropy(real p, real n)
{
    return -p * log2(p) - n * log2(n);
}

real entropy(real p)
{
    // TODO: could be more accurate by using log1p(-p) for log(1-p)
    return (p == 0 ? 0 : -p * log2(p)) -
           (p == 1 ? 0 : (1 - p) * log2(1 - p));
}

/+real entropy(alias pred, Range)(Range r)
    if (isInputRange!Range)
{
    size_t pos = 0;
    size_t total = 0;
    foreach (e; r)
    {
        if (pred(e))
            ++pos;
        ++total;
    }
    return entropy(cast(real)pos / total);
}+/

real binaryAttributeInfoGain(alias classPred, alias attributePred, Range)(Range r)
    if (isForwardRange!Range)
{
    size_t nT1 = 0;
    size_t nT2 = 0;
    size_t nF1 = 0;
    size_t nF2 = 0;
    foreach (e; r)
    {
        if (attributePred(e))
        {
            if (classPred(e))
                nT1++;
            else
                nT2++;
        }
        else
        {
            if (classPred(e))
                nF1++;
            else 
                nF2++;
        }
    }
    size_t nT = nT1 + nT2;
    size_t nF = nF1 + nF2;
    size_t n1 = nT1 + nF1;
    size_t n2 = nT2 + nF2;
    size_t N = nT + nF;
    if (N == 0)
        return 0.0;
    return entropy(cast(real)n1 / N) -
           (nT == 0 ? 0 : entropy(cast(real)nT1 / nT) * nT / N) -
           (nF == 0 ? 0 : entropy(cast(real)nF1 / nF) * nF / N);
}

__gshared uint[] nT1Memory;
__gshared uint[] nT2Memory;

struct Info(AttributeType)
{
    AttributeType attribute;
    double infoGainProxy;
    double t1ratio;
    double f1ratio;
    size_t nT1, nF1, nT2, nF2, nT, nF, n1, n2;
}

// Discovered through grid search
size_t ttermThres = 20;
size_t ftermThres = 10;

auto attributeInfoGain(alias classPred, alias attributePred, TrainingRange, AttributeRange)
                      (TrainingRange trainingSet, AttributeRange attributes)
{
    size_t nAttribs = attributes.save.walkLength;
    if (nT1Memory.length < nAttribs)
        nT1Memory = new uint[nAttribs];
    if (nT2Memory.length < nAttribs)
        nT2Memory = new uint[nAttribs];
    uint[] nT1 = nT1Memory;
    uint[] nT2 = nT2Memory;
    nT1[] = 0;
    nT2[] = 0;
    size_t N = 0;
    size_t n1 = 0;
    foreach (sample; trainingSet)
    {
        auto isClass1 = classPred(sample);
        /+size_t i = 0;
        foreach (attribute; attributes.save)
        {
            if (attributePred(sample, attribute))
            {
                if (isClass1)
                    nT1[i]++;
                else
                    nT2[i]++;
            }
            ++i;
        }+/
        // HACK
        foreach (t; sample.terms)
            if (isClass1)
                nT1[t]++;
            else
                nT2[t]++;


        if (isClass1)
            n1++;
        ++N;
    }
    size_t n2 = N - n1;

    auto best = Info!(ElementType!AttributeRange)(ushort.max, -9999999.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0);
    size_t i = 0;
    foreach (attribute; attributes)
    {
        size_t nT = nT1[i] + nT2[i];
        size_t cnF1 = n1 - nT1[i];
        size_t cnF2 = n2 - nT2[i];
        size_t nF = cnF1 + cnF2;
        if (nT > ttermThres && nF > ftermThres)
        {
            /+real infoGain = entropy(cast(real)n1 / N) -
                   (nT == 0 ? 0 : entropy(cast(real)nT1[i] / nT) * nT / N) -
                   (nF == 0 ? 0 : entropy(cast(real)cnF1 / nF) * nF / N);+/
            real infoGainProxy = -entropy(cast(real)nT1[i] / nT) * nT -
                                  entropy(cast(real)cnF1 / nF) * nF;
            if (infoGainProxy > best.infoGainProxy)
                best = Info!(ElementType!AttributeRange)(
                             attribute,
                             cast(double)infoGainProxy,
                             cast(double)nT1[i] / nT,
                             cast(double)cnF1 / nF,
                             nT1[i], cnF1, nT2[i], cnF2,
                             nT, nF, n1, n2 );
        }
        ++i;
    }
    return best;
}