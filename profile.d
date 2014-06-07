module stdex.profile;

// Only enable this code in a profile build.
version(profile):

import stdex.ucontext;
import core.sys.posix.time;
import core.time;
import core.stdc.time;
import core.sys.posix.config;
import core.sys.posix.signal;
import std.c.stdio;
import std.c.stdlib;
import std.c.string;
import std.range;
import std.stdio;
import std.demangle;
import std.exception;

extern(C) struct timespec
{
    time_t tv_sec;
    c_long tv_nsec;
}

extern(C) struct itimerval
{
    timespec it_interval;
    timespec it_value;
}

enum ITIMER_REAL = 0;
enum ITIMER_VIRTUAL = 1;
enum ITIMER_PROF = 2;

extern(C) int setitimer(int, const(itimerval)*, itimerval*);

shared static this()
{
    g_samples = (cast(Record*)malloc((1<<20) * Record.sizeof))[0..1<<20];

    // SIGPROF setup
    sigaction_t sa;
    memset(&sa, 0, sa.sizeof);
    sa.sa_sigaction = &handler;
    sa.sa_flags = SA_RESTART;
    memset(&sa.sa_mask, 0xff, sa.sa_mask.sizeof);
    sigaction(SIGPROF, &sa, null);

    // setup profile timer (1000Hz)
    itimerval it;
    memset(&it, 0, it.sizeof);
    it.it_interval.tv_nsec = 1000;
    it.it_value = it.it_interval;
    setitimer(ITIMER_PROF, &it, null);

    // setup SIGINT handler
    sigaction_t sa2;
    memset(&sa2, 0, sa2.sizeof);
    sa2.sa_sigaction = &sigintHandler;
    //sa2.sa_flags = SA_RESTART;
    memset(&sa2.sa_mask, 0xff, sa2.sa_mask.sizeof);
    sigaction(SIGINT, &sa2, null);
}

shared static ~this()
{
    flush(File("profile_output", "w"));
}

extern(C) void sigintHandler(int sig, siginfo_t* siginfo, void* context)
{
    flush(File("profile_output", "w"));
    exit(1);
}

extern(C) int backtrace(void** buffer, int size);
extern(C) char** backtrace_symbols(void** buffer, int size);

extern(C) void handler(int sig, siginfo_t* siginfo, void* context)
{
    ucontext* ucon = cast(ucontext*)context;
    void* frame = cast(void*)ucon.uc_mcontext.__ss.__rip;
    enum skip = 3;
    void*[19+skip] temp;
    g_samples[g_written].pc[0] = frame;
    g_samples[g_written].n = backtrace(temp.ptr, cast(int)temp.length) - skip + 1;
    g_samples[g_written].pc[1..$] = temp[skip..$];
    g_written++;
}

struct Record
{
    void*[20] pc;
    size_t n;
}

__gshared Record[] g_samples;
__gshared size_t g_written = 0;
__gshared size_t g_read = 0;

Record[] newSamples()
{
    return g_samples[g_read..g_written];
}

void flush(File file)
{
    writeln("Flushing profile samples...");
    auto psamples = newSamples();
    foreach (ref rec; psamples)
    {
        foreach (pc; rec.pc[0..1])//rec.n])
        {
            char** sym = backtrace_symbols(&pc, 1);
            string ss = cast(string)sym[0][0..strlen(sym[0])];
            file.writeln(demangle(ss.split.array[3]));
            free(sym);
        }
    }
}