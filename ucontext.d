module stdex.ucontext;

alias uint sigset_t;

// For __pc
enum FP_PREC_24B = 0;
enum FP_PREC_53B = 2;
enum FP_PREC_64B = 3;

// For __rc
enum FP_RND_NEAR = 0;
enum FP_RND_DOWN = 1;
enum FP_RND_UP = 2;
enum FP_CHOP = 3;

extern(C) struct fp_control
{
    private ushort bits;

    @property
    {
        ushort __invalid() { return bits & 1; }
        ushort __denorm() { return (bits >> 1) & 1; }
        ushort __zdiv() { return (bits >> 2) & 1; }
        ushort __ovrfl() { return (bits >> 3) & 1; }
        ushort __undfl() { return (bits >> 4) & 1; }
        ushort __precis() { return (bits >> 5) & 1; }
        ushort __pc() { return (bits >> 8) & 2; }
        ushort __rc() { return (bits >> 10) & 2; }
    }
    /+  ushort  __invalid   :1,
                __denorm    :1,
                __zdiv      :1,
                __ovrfl     :1,
                __undfl     :1,
                __precis    :1,
                        :2,
                __pc        :2,
                __rc        :2,
                /*inf*/     :1,
                        :3;     +/
}

extern(C) struct fp_status
{
    private ushort bits;

    @property
    {
        ushort __invalid() { return bits & 1; }
        ushort __denorm() { return (bits >> 1) & 1; }
        ushort __zdiv() { return (bits >> 2) & 1; }
        ushort __ovrfl() { return (bits >> 3) & 1; }
        ushort __undfl() { return (bits >> 4) & 1; }
        ushort __precis() { return (bits >> 5) & 1; }
        ushort __stkflt() { return (bits >> 6) & 1; }
        ushort __errsumm() { return (bits >> 7) & 1; }
        ushort __c0() { return (bits >> 8) & 1; }
        ushort __c1() { return (bits >> 9) & 1; }
        ushort __c2() { return (bits >> 10) & 1; }
        ushort __tos() { return (bits >> 11) & 3; }
        ushort __c3() { return (bits >> 14) & 1; }
        ushort __busy() { return (bits >> 15) & 1; }
    }
    /+  ushort  __invalid   :1,
                __denorm    :1,
                __zdiv      :1,
                __ovrfl     :1,
                __undfl     :1,
                __precis    :1,
                __stkflt    :1,
                __errsumm   :1,
                __c0        :1,
                __c1        :1,
                __c2        :1,
                __tos       :3,
                __c3        :1,
                __busy      :1; +/
}

struct __darwin_mmst_reg
{
    char __mmst_reg[10];
    char __mmst_rsrv[6];
}

struct __darwin_xmm_reg
{
    char __xmm_reg[16];
}

enum FP_STATE_BYTES = 512; /* number of chars worth of data from fpu_fcw */

version (X86)
{   
    extern(C) struct thread_state
    {
        uint __eax;
        uint __ebx;
        uint __ecx;
        uint __edx;
        uint __edi;
        uint __esi;
        uint __ebp;
        uint __esp;
        uint __ss;
        uint __eflags;
        uint __eip;
        uint __cs;
        uint __ds;
        uint __es;
        uint __fs;
        uint __gs;
    }

    extern(C) struct float_state
    {
        int __fpu_reserved[2];
        fp_control __fpu_fcw;      /* x87 FPU control word */
        fp_status __fpu_fsw;      /* x87 FPU status word */
        ubyte __fpu_ftw;      /* x87 FPU tag word */
        ubyte __fpu_rsrv1;        /* reserved */ 
        ushort __fpu_fop;      /* x87 FPU Opcode */
        uint __fpu_ip;       /* x87 FPU Instruction Pointer offset */
        ushort __fpu_cs;       /* x87 FPU Instruction Pointer Selector */
        ushort __fpu_rsrv2;        /* reserved */
        uint __fpu_dp;       /* x87 FPU Instruction Operand(Data) Pointer offset */
        ushort __fpu_ds;       /* x87 FPU Instruction Operand(Data) Pointer Selector */
        ushort __fpu_rsrv3;        /* reserved */
        uint __fpu_mxcsr;        /* MXCSR Register state */
        uint __fpu_mxcsrmask;    /* MXCSR mask */
        __darwin_mmst_reg __fpu_stmm0;        /* ST0/MM0   */
        __darwin_mmst_reg __fpu_stmm1;        /* ST1/MM1  */
        __darwin_mmst_reg __fpu_stmm2;        /* ST2/MM2  */
        __darwin_mmst_reg __fpu_stmm3;        /* ST3/MM3  */
        __darwin_mmst_reg __fpu_stmm4;        /* ST4/MM4  */
        __darwin_mmst_reg __fpu_stmm5;        /* ST5/MM5  */
        __darwin_mmst_reg __fpu_stmm6;        /* ST6/MM6  */
        __darwin_mmst_reg __fpu_stmm7;        /* ST7/MM7  */
        __darwin_xmm_reg __fpu_xmm0;     /* XMM 0  */
        __darwin_xmm_reg __fpu_xmm1;     /* XMM 1  */
        __darwin_xmm_reg __fpu_xmm2;     /* XMM 2  */
        __darwin_xmm_reg __fpu_xmm3;     /* XMM 3  */
        __darwin_xmm_reg __fpu_xmm4;     /* XMM 4  */
        __darwin_xmm_reg __fpu_xmm5;     /* XMM 5  */
        __darwin_xmm_reg __fpu_xmm6;     /* XMM 6  */
        __darwin_xmm_reg __fpu_xmm7;     /* XMM 7  */
        char __fpu_rsrv4[14*16]; /* reserved */
        int __fpu_reserved1;
    }

    extern(C) struct exception_state
    {
        uint __trapno;
        uint __err;
        uint __faultvaddr;
    }

    extern(C) struct debug_state
    {
        uint __dr0;
        uint __dr1;
        uint __dr2;
        uint __dr3;
        uint __dr4;
        uint __dr5;
        uint __dr6;
        uint __dr7;
    }
}
else version (X86_64)
{
    extern(C) struct thread_state
    {
        ulong __rax;
        ulong __rbx;
        ulong __rcx;
        ulong __rdx;
        ulong __rdi;
        ulong __rsi;
        ulong __rbp;
        ulong __rsp;
        ulong __r8;
        ulong __r9;
        ulong __r10;
        ulong __r11;
        ulong __r12;
        ulong __r13;
        ulong __r14;
        ulong __r15;
        ulong __rip;
        ulong __rflags;
        ulong __cs;
        ulong __fs;
        ulong __gs;
    }

    extern(C) struct float_state
    {
        int __fpu_reserved[2];
        fp_control __fpu_fcw;      /* x87 FPU control word */
        fp_status __fpu_fsw;      /* x87 FPU status word */
        ubyte __fpu_ftw;      /* x87 FPU tag word */
        ubyte __fpu_rsrv1;        /* reserved */ 
        ushort __fpu_fop;      /* x87 FPU Opcode */

        /* x87 FPU Instruction Pointer */
        uint __fpu_ip;       /* offset */
        ushort __fpu_cs;       /* Selector */

        ushort __fpu_rsrv2;        /* reserved */

        /* x87 FPU Instruction Operand(Data) Pointer */
        uint __fpu_dp;       /* offset */
        ushort __fpu_ds;       /* Selector */

        ushort __fpu_rsrv3;        /* reserved */
        uint __fpu_mxcsr;        /* MXCSR Register state */
        uint __fpu_mxcsrmask;    /* MXCSR mask */
        __darwin_mmst_reg __fpu_stmm0;        /* ST0/MM0   */
        __darwin_mmst_reg __fpu_stmm1;        /* ST1/MM1  */
        __darwin_mmst_reg __fpu_stmm2;        /* ST2/MM2  */
        __darwin_mmst_reg __fpu_stmm3;        /* ST3/MM3  */
        __darwin_mmst_reg __fpu_stmm4;        /* ST4/MM4  */
        __darwin_mmst_reg __fpu_stmm5;        /* ST5/MM5  */
        __darwin_mmst_reg __fpu_stmm6;        /* ST6/MM6  */
        __darwin_mmst_reg __fpu_stmm7;        /* ST7/MM7  */
        __darwin_xmm_reg __fpu_xmm0;     /* XMM 0  */
        __darwin_xmm_reg __fpu_xmm1;     /* XMM 1  */
        __darwin_xmm_reg __fpu_xmm2;     /* XMM 2  */
        __darwin_xmm_reg __fpu_xmm3;     /* XMM 3  */
        __darwin_xmm_reg __fpu_xmm4;     /* XMM 4  */
        __darwin_xmm_reg __fpu_xmm5;     /* XMM 5  */
        __darwin_xmm_reg __fpu_xmm6;     /* XMM 6  */
        __darwin_xmm_reg __fpu_xmm7;     /* XMM 7  */
        __darwin_xmm_reg __fpu_xmm8;     /* XMM 8  */
        __darwin_xmm_reg __fpu_xmm9;     /* XMM 9  */
        __darwin_xmm_reg __fpu_xmm10;        /* XMM 10  */
        __darwin_xmm_reg __fpu_xmm11;        /* XMM 11 */
        __darwin_xmm_reg __fpu_xmm12;        /* XMM 12  */
        __darwin_xmm_reg __fpu_xmm13;        /* XMM 13  */
        __darwin_xmm_reg __fpu_xmm14;        /* XMM 14  */
        __darwin_xmm_reg __fpu_xmm15;        /* XMM 15  */
        char __fpu_rsrv4[6*16];  /* reserved */
        int __fpu_reserved1;
    }

    extern(C) struct exception_state
    {
        uint __trapno;
        uint __err;
        ulong __faultvaddr;
    }

    extern(C) struct debug_state
    {
        ulong __dr0;
        ulong __dr1;
        ulong __dr2;
        ulong __dr3;
        ulong __dr4;
        ulong __dr5;
        ulong __dr6;
        ulong __dr7;
    }
}
else
{
    static assert("Platform not supported.");
}

extern(C) struct mcontext
{
    exception_state __es;
    thread_state __ss;
    float_state __fs;
}

extern(C) struct sigaltstack
{
    void *ss_sp;        /* signal stack base */
    size_t ss_size;     /* signal stack length */
    int ss_flags;       /* SA_DISABLE and/or SA_ONSTACK */
}

extern(C) struct ucontext
{
    int uc_onstack;
    sigset_t uc_sigmask; /* signal mask used by this context */
    sigaltstack uc_stack;   /* stack used by this context */
    ucontext* uc_link;   /* pointer to resuming context */
    size_t uc_mcsize;  /* size of the machine context passed in */
    mcontext* uc_mcontext;   /* pointer to machine specific context */
}
