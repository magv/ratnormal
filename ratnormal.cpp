static const char usagetext[] = R"(
Ss{NAME}
    Nm{ratnormal} -- compute rational normal form of a sum of rational
    expressions.

Ss{SYNOPSYS}
    Nm{ratnormal} [Fl{-hC}] [Fl{-j} Ar{threads}] Ar{input} [Ar{output}]

Ss{DESCRIPTION}

    Nm{ratnormal} reads an input file, which must contain a sum of
    rational expressions (products of integer powers of polynomials with
    integer coefficients), and brings the sum under a common denominator,
    cancelling the fractions if needed, and preserving the factorization
    of the denominator.
    
    Nm{ratnormal} uses the FLINT library.

Ss{EXAMPLES}

    $ echo '2/3/(x + 2*y^2)^2 - 1/(x + 6*x*y^2)' | ratnormal -
    (
     1/3 *
     1/(6*x*y^2+x) *
     1/(x+2*y^2)^2 *
     (-3*x^2+2*x-12*y^4)
    )

Ss{OPTIONS}
    Fl{-j} Ar{threads} Use this many worker threads, if possible (default: 1).
    Fl{-h}         Show this help message.
    Fl{-C}         Force colored output even if stderr is not a tty.

Ss{ARGUMENTS}
    Ar{input}      Read the input expression from this file, with "-"
               meaning the standard input.
    Ar{output}     Write the result into this file, with "-" meaning
               the standard output (which is the default).

Ss{AUTHORS}
    Vitaly Magerya <vitaly.magerya@tx97.net>
)";

#include <assert.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <fstream>
#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <iostream>
#include <printf.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

using namespace std;

static bool COLORS = !!isatty(STDERR_FILENO);

/* Misc
 */

#define SWAP(type, a, b) { type t = a; a = b; b = t; }

static double
timestamp(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

/* Logging
 */

static int log_depth = 0;
static double log_first_timestamp = timestamp();
static double log_last_timestamp = log_first_timestamp;

typedef struct {
    const char *name;
    double time;
} _log_func_info;

#define logd(fmt, ...) {\
    double NOW = timestamp(); \
    fprintf(stderr, "[%.4fs +%.4fs]%*s " fmt "\n", NOW - log_first_timestamp, NOW - log_last_timestamp, log_depth*2, "", ## __VA_ARGS__); \
    log_last_timestamp = NOW; \
}

void
_log_func_start(const _log_func_info *i)
{
    logd("> %s()", i->name);
    log_depth++;
}

void
_log_func_end(const _log_func_info *i)
{
    log_depth--;
    logd("< %s(%.4fs)", i->name, NOW - i->time);
}

#define LOGME \
    __attribute__((cleanup(_log_func_end))) _log_func_info _log_func_info = {__func__, timestamp()}; \
    _log_func_start(&_log_func_info);

/* Rational context
 */

// struct rat_ctx {
fmpz_mpoly_ctx_t ctx;
map<string, int> variable2index;
vector<string> variables;
const char **variable_names;
#define nvariables (variables.size())
// }

void
rat_ctx_init()
{
    fmpz_mpoly_ctx_init(ctx, nvariables, ORD_LEX);
}

void
rat_ctx_clear()
{
    fmpz_mpoly_ctx_clear(ctx);
    for (ulong i = 0; i < nvariables; i++) {
        flint_free((void*)variable_names[i]);
        variable_names[i] = NULL;
    }
}

/* Rational
 */

typedef struct {
    fmpq_t numfactor;
    fmpz_mpoly_struct *factors;
    int *powers;
    ulong num;
    ulong alloc;
} rat_struct;

typedef rat_struct rat_t[1];

void
rat_init(rat_t rat)
{
    fmpq_init(rat->numfactor);
    rat->factors = NULL;
    rat->powers = NULL;
    rat->num = 0;
    rat->alloc = 0;
}

void
rat_one(rat_t rat)
{
    fmpq_one(rat->numfactor);
    for (ulong i = 0; i < rat->num; i++) {
        fmpz_mpoly_clear(&rat->factors[i], ctx);
    }
    rat->num = 0;
}

void
rat_swap(rat_t a, rat_t b)
{
   rat_struct t = *a;
   *a = *b;
   *b = t;
}

void
rat_clear(rat_t rat)
{
    fmpq_clear(rat->numfactor);
    for (ulong i = 0; i < rat->num; i++) {
        fmpz_mpoly_clear(&rat->factors[i], ctx);
    }
    flint_free(rat->factors);
    flint_free(rat->powers);
    // if paranoid then {
    rat->factors = NULL;
    rat->powers = NULL;
    rat->num = 0;
    rat->alloc = 0;
    // }
}

void
rat_fit_length(rat_t rat, ulong len)
{
    if (rat->alloc < len) {
        len = FLINT_MAX(len, rat->alloc + rat->alloc/2);
        rat->factors = (fmpz_mpoly_struct*)flint_realloc(rat->factors, len*sizeof(rat->factors[0]));
        rat->powers = (int*)flint_realloc(rat->powers, len*sizeof(rat->powers[0]));
        rat->alloc = len;
    }
}

void
rat_mul_fmpq(rat_t rat, const fmpq_t x, const int power)
{
    assert(!fmpq_is_zero(x));
    fmpq_t xpow;
    fmpq_init(xpow);
    fmpq_pow_si(xpow, x, power);
    fmpq_mul(rat->numfactor, rat->numfactor, xpow);
    fmpq_clear(xpow);
}

void
rat_mul_fmpz(rat_t rat, const fmpz_t x, const int power)
{
    assert(!fmpz_is_zero(x));
    if (power == 0) {
        return;
    } else {
        fmpz_t xpow;
        fmpz_init(xpow);
        if (power > 0) {
            fmpz_pow_ui(xpow, x, power);
            fmpq_mul_fmpz(rat->numfactor, rat->numfactor, xpow);
        } else {
            fmpz_pow_ui(xpow, x, -power);
            fmpq_div_fmpz(rat->numfactor, rat->numfactor, xpow);
        }
        fmpz_clear(xpow);
    }
}

void
rat_mul_fmpz_mpoly(rat_t rat, const fmpz_mpoly_t poly, int power)
{
    if (power == 0) {
        return;
    } else if (fmpz_mpoly_is_fmpz(poly, ctx)) {
        fmpz_t k;
        fmpz_mpoly_get_fmpz(k, poly, ctx);
        rat_mul_fmpz(rat, k, power);
    } else {
        rat_fit_length(rat, rat->num + 1);
        fmpz_mpoly_init(&rat->factors[rat->num], ctx);
        fmpz_mpoly_set(&rat->factors[rat->num], poly, ctx);
        rat->powers[rat->num] = power;
        rat->num++;
    }
}

void
rat_mul_fmpz_mpoly_setx(rat_t rat, fmpz_mpoly_t poly, int power)
{
    if (power == 0) {
        return;
    } else if (fmpz_mpoly_is_fmpz(poly, ctx)) {
        fmpz_t k;
        fmpz_init(k);
        fmpz_mpoly_get_fmpz(k, poly, ctx);
        rat_mul_fmpz(rat, k, power);
        fmpz_mpoly_zero(poly, ctx);
    } else {
        rat_fit_length(rat, rat->num + 1);
        fmpz_mpoly_init(&rat->factors[rat->num], ctx);
        fmpz_mpoly_swap(&rat->factors[rat->num], poly, ctx);
        rat->powers[rat->num] = power;
        rat->num++;
    }
}

void
rat_fprint(FILE *f, const rat_t rat)
{
    if ((rat->num == 0) || !fmpq_is_one(rat->numfactor)) {
        fprintf(f, " ");
        fmpq_fprint(f, rat->numfactor);
        if (0 != rat->num) {
            fprintf(f, " *\n");
        }
    }
    for (ulong i = 0; i < rat->num; i++) {
        int power = rat->powers[i];
        if (power >= 0) {
            fprintf(f, " (");
        } else {
            fprintf(f, " 1/(");
            power = -power;
        }
        fmpz_mpoly_fprint_pretty(f, &rat->factors[i], variable_names, ctx);
        if (power != 1) {
            fprintf(f, ")^%d", power);
        } else {
            fprintf(f, ")");
        }
        if (i != rat->num-1) {
            fprintf(f, " *\n");
        }
    }
}

void
rat_fprint_short(FILE *f, const rat_t rat)
{
    fprintf(f, "R{");
    for (ulong i = 0; i < rat->num; i++) {
        slong nterms = fmpz_mpoly_length(&rat->factors[i], ctx);
        int power = rat->powers[i];
        if (i == 0) {
            if (power < 0) {
                fprintf(f, "1/");
                power = -power;
            }
        } else {
            if (power < 0) {
                fprintf(f, "/");
                power = -power;
            } else {
                fprintf(f, " ");
            }
        }
        if (power == 1) {
            fprintf(f, "%ldt", nterms);
        } else {
            fprintf(f, "%ldt^%d", nterms, power);
        }
    }
    fprintf(f, "}");
}

void
rat_swapoff_fmpz_factor(rat_t rat, ulong i)
{
    assert(i < rat->num);
    assert(fmpz_mpoly_is_fmpz(&rat->factors[i], ctx));
    rat_mul_fmpz_mpoly_setx(rat, &rat->factors[i], rat->powers[i]);
    if (i < rat->num-1) {
        fmpz_mpoly_swap(&rat->factors[i], &rat->factors[rat->num-1], ctx);
        SWAP(int, rat->powers[i], rat->powers[rat->num-1]);
    }
    fmpz_mpoly_clear(&rat->factors[rat->num-1], ctx);
    rat->num--;
}

slong
rat_max_length(const rat_t rat)
{
    slong maxlen = 0;
    for (ulong i = 0; i < rat->num; i++) {
        maxlen = FLINT_MAX(maxlen, fmpz_mpoly_length(&rat->factors[i], ctx));
    }
    return maxlen;
}

void
rat_sort(rat_t rat)
{
    for (ulong i = 0; i < rat->num; i++) {
        slong nterms1 = fmpz_mpoly_length(&rat->factors[i], ctx);
        for (ulong j = i+1; j < rat->num; j++) {
            slong nterms2 = fmpz_mpoly_length(&rat->factors[j], ctx);
            if (nterms2 > nterms1) {
                fmpz_mpoly_swap(&rat->factors[i], &rat->factors[j], ctx);
                SWAP(int, rat->powers[i], rat->powers[j]);
                nterms1 = nterms2;
            }
        }
    }
}

void
rat_reverse(rat_t rat)
{
    for (ulong i = 0; i < rat->num/2; i++) {
        ulong j = rat->num - 1 - i;
        assert(i < j);
        fmpz_mpoly_swap(&rat->factors[i], &rat->factors[j], ctx);
        SWAP(int, rat->powers[i], rat->powers[j]);
    }
}

void
rat_cofactorize(rat_t rat)
{
    LOGME;
    logd("Cofactorizing %r", rat);
    rat_sort(rat);
    fmpz_mpoly_t gcd, aprime, bprime;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(aprime, ctx);
    fmpz_mpoly_init(bprime, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        #define A (&rat->factors[i])
        #define Apower rat->powers[i]
        //if (Apower >= 0) continue;
        for (ulong j = i + 1; j < rat->num; j++) {
            #define B (&rat->factors[j])
            #define Bpower rat->powers[j]
            //if (Bpower >= 0) continue;
            fmpz_mpoly_gcd(gcd, A, B, ctx);
            if (fmpz_mpoly_is_one(gcd, ctx)) {
                continue;
            }
            logd("Split factor: %p^%d from %p^%d %p^%d", gcd, Apower + Bpower, A, Apower, B, Bpower);
            if (fmpz_mpoly_divides(aprime, A, gcd, ctx) != 1) { assert(0); }
            if (fmpz_mpoly_divides(bprime, B, gcd, ctx) != 1) { assert(0); }
            fmpz_mpoly_swap(A, aprime, ctx);
            fmpz_mpoly_swap(B, bprime, ctx);
            rat_mul_fmpz_mpoly_setx(rat, gcd, Apower + Bpower);
            if (fmpz_mpoly_is_fmpz(B, ctx)) {
                rat_swapoff_fmpz_factor(rat, j);
                j--;
            }
            if (fmpz_mpoly_is_fmpz(A, ctx)) {
                rat_swapoff_fmpz_factor(rat, i);
                i--;
                break;
            }
            #undef B
            #undef Bpower
        }
        #undef A
        #undef Apower
    }
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(aprime, ctx);
    fmpz_mpoly_clear(bprime, ctx);
    rat_sort(rat);
}

void
rat_expand_numerator(fmpz_mpoly_t poly, const rat_t rat)
{
    fmpz_mpoly_t tmp;
    fmpz_mpoly_init(tmp, ctx);
    fmpz_mpoly_one(poly, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        if (rat->powers[i] <= 0) continue;
        if (rat->powers[i] > 1) {
            fmpz_mpoly_pow_ui(tmp, &rat->factors[i], rat->powers[i], ctx);
            fmpz_mpoly_mul(poly, poly, tmp, ctx);
        } else {
            fmpz_mpoly_mul(poly, poly, &rat->factors[i], ctx);
        }
    }
    fmpz_mpoly_clear(tmp, ctx);
    if (!fmpz_is_one(fmpq_numref(rat->numfactor))) {
        fmpz_mpoly_scalar_mul_fmpz(poly, poly, fmpq_numref(rat->numfactor), ctx);
    }
}

void
rat_expand_denominator(fmpz_mpoly_t poly, const rat_t rat)
{
    fmpz_mpoly_t tmp;
    fmpz_mpoly_init(tmp, ctx);
    fmpz_mpoly_one(poly, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        if (rat->powers[i] >= 0) continue;
        if (rat->powers[i] < -1) {
            fmpz_mpoly_pow_ui(tmp, &rat->factors[i], -rat->powers[i], ctx);
            fmpz_mpoly_mul(poly, poly, tmp, ctx);
        } else {
            fmpz_mpoly_mul(poly, poly, &rat->factors[i], ctx);
        }
    }
    fmpz_mpoly_clear(tmp, ctx);
    if (!fmpz_is_one(fmpq_denref(rat->numfactor))) {
        fmpz_mpoly_scalar_mul_fmpz(poly, poly, fmpq_denref(rat->numfactor), ctx);
    }
}

void
rat_add_setx(rat_t res, rat_t rat1, rat_t rat2)
{
    LOGME;
    logd("Adding %r and %r", rat1, rat2);
    rat_one(res);
    // First, take out common factors out of the denominators.
    fmpz_mpoly_t gcd, aprime, bprime;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(aprime, ctx);
    fmpz_mpoly_init(bprime, ctx);
    #define A (&rat1->factors[i])
    #define Apower rat1->powers[i]
    #define B (&rat2->factors[j])
    #define Bpower rat2->powers[j]
    for (ulong i = 0; i < rat1->num; i++) {
        if (Apower >= 0) continue;
        for (ulong j = 0; j < rat2->num; j++) {
            if (Bpower >= 0) continue;
            fmpz_mpoly_gcd_cofactors(gcd, aprime, bprime, A, B, ctx);
            if (fmpz_mpoly_is_one(gcd, ctx)) {
                continue;
            }
            int power = FLINT_MAX(Apower, Bpower);
            //logd("Common factor: %p^%d", gcd, power);
            fmpz_mpoly_swap(A, aprime, ctx);
            fmpz_mpoly_swap(B, bprime, ctx);
            rat_mul_fmpz_mpoly(rat1, gcd, Apower - power);
            rat_mul_fmpz_mpoly(rat2, gcd, Bpower - power);
            rat_mul_fmpz_mpoly_setx(res, gcd, power);
            if (fmpz_mpoly_is_fmpz(B, ctx)) {
                rat_swapoff_fmpz_factor(rat2, j);
                j--;
            }
            if (fmpz_mpoly_is_fmpz(A, ctx)) {
                rat_swapoff_fmpz_factor(rat1, i);
                i--;
                break;
            }
        }
    }
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(aprime, ctx);
    fmpz_mpoly_clear(bprime, ctx);
    // Don't forget about numeric factors.
    fmpz_t g;
    fmpz_init(g);
    fmpz_gcd(g, fmpq_denref(rat1->numfactor), fmpq_denref(rat2->numfactor));
    fmpq_mul_fmpz(rat1->numfactor, rat1->numfactor, g);
    fmpq_mul_fmpz(rat2->numfactor, rat2->numfactor, g);
    rat_mul_fmpz(res, g, -1);
    fmpz_clear(g);
    logd("Common factor: %r", res);
    // Then, expand the combined numerator.
    fmpz_mpoly_t anum, aden, bnum, bden;
    fmpz_mpoly_init(anum, ctx);
    fmpz_mpoly_init(bnum, ctx);
    fmpz_mpoly_init(aden, ctx);
    fmpz_mpoly_init(bden, ctx);
    rat_expand_numerator(anum, rat1);
    rat_expand_denominator(aden, rat1);
    rat_expand_numerator(bnum, rat2);
    rat_expand_denominator(bden, rat2);
    logd("Computing %p*%p + %p*%p", anum, bden, bnum, aden);
    fmpz_mpoly_mul(anum, anum, bden, ctx);
    fmpz_mpoly_mul(bnum, bnum, aden, ctx);
    fmpz_mpoly_add(anum, anum, bnum, ctx);
    rat_mul_fmpz_mpoly_setx(res, anum, 1);
    fmpz_mpoly_clear(anum, ctx);
    fmpz_mpoly_clear(aden, ctx);
    fmpz_mpoly_clear(bnum, ctx);
    fmpz_mpoly_clear(bden, ctx);
    // Finally, multiply by the combined denominator.
    for (ulong i = 0; i < rat1->num; i++) {
        if (Apower >= 0) continue;
        rat_mul_fmpz_mpoly_setx(res, A, Apower);
    }
    for (ulong j = 0; j < rat2->num; j++) {
        if (Bpower >= 0) continue;
        rat_mul_fmpz_mpoly_setx(res, B, Bpower);
    }
    rat_mul_fmpz(res, fmpq_denref(rat1->numfactor), -1);
    rat_mul_fmpz(res, fmpq_denref(rat2->numfactor), -1);
    #undef A
    #undef Apower
    #undef B
    #undef Bpower
    rat_cofactorize(res);
    logd("Result is %r", res);
}

/* Rational sum
 */

typedef struct {
    rat_struct *terms;
    ulong num;
    ulong alloc;
} ratsum_struct;

typedef ratsum_struct ratsum_t[1];

void
ratsum_init(ratsum_t sum)
{
    sum->terms = NULL;
    sum->num = 0;
    sum->alloc = 0;
}

void
ratsum_zero(ratsum_t sum)
{
    for (ulong i = 0; i < sum->num; i++) {
        rat_clear(&sum->terms[i]);
    }
    sum->num = 0;
}

void
ratsum_fit_length(ratsum_t sum, ulong len)
{
    if (sum->alloc < len) {
        len = FLINT_MAX(len, sum->alloc + sum->alloc/2);
        sum->terms = (rat_struct*)flint_realloc(sum->terms, len*sizeof(sum->terms[0]));
        sum->alloc = len;
    }
}

void
ratsum_clear(ratsum_t sum)
{
    for (ulong i = 0; i < sum->num; i++) {
        rat_clear(&sum->terms[i]);
    }
    flint_free(sum->terms);
    sum->terms = NULL;
    sum->num = 0;
    sum->alloc = 0;
}

void
ratsum_fprint(FILE *f, const ratsum_t sum)
{
    fprintf(f, "(\n ");
    for (ulong i = 0; i < sum->num; i++) {
        if (i != 0) fprintf(f, " +\n ");
        rat_fprint(f, &sum->terms[i]);
    }
    fprintf(f, "\n)");
}

void
ratsum_add_rat_setx(ratsum_t sum, rat_t rat)
{
    ratsum_fit_length(sum, sum->num+1);
    rat_init(&sum->terms[sum->num]);
    rat_swap(&sum->terms[sum->num], rat);
    sum->num++;
}

void
ratsum_sum_setx(rat_t rat, ratsum_t sum)
{
    LOGME;
    for (ulong n = 0; n < sum->num-1; n++) {
        // Find two shortest rationals, add them.
        slong length1 = LONG_MAX, idx1 = -1;
        slong length2 = LONG_MAX, idx2 = -1;
        for (ulong i = 0; i < sum->num - n; i++) {
            slong length = rat_max_length(&sum->terms[i]);
            if (length < length1) {
                length2 = length1; idx2 = idx1;
                length1 = length; idx1 = i;
                continue;
            }
            if (length < length2) {
                length2 = length; idx2 = i;
                continue;
            }
        }
        assert((idx1 >= 0) && (idx2 >= 0) && (idx1 != idx2));
        logd("Adding pair %ld of %ld", n+1, sum->num-1);
        rat_add_setx(rat, &sum->terms[idx1], &sum->terms[idx2]);
        rat_swap(&sum->terms[idx1], rat);
        rat_swap(&sum->terms[idx2], &sum->terms[sum->num - n - 1]);
    }
    rat_swap(&sum->terms[0], rat);
}

/* GiNaC conversion
 */

using namespace GiNaC;

template <typename F> void
factor_iter(const ex &e, F yield)
{
    if (is_a<mul>(e)) {
        for (const auto &f : e) {
            if (is_a<power>(f)) {
                yield(f.op(0), ex_to<numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (is_a<power>(e)) {
            yield(e.op(0), ex_to<numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

template <typename F> void
term_iter(const ex &e, F yield)
{
    if (is_a<add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

void
fmpz_of_ginac(fmpz_t x, const GiNaC::numeric &num)
{
    assert(num.is_integer());
    assert((LONG_MIN <= num) && (num <= LONG_MAX));
    fmpz_set_si(x, num.to_long());
}

void
fmpq_of_ginac(fmpq_t x, const GiNaC::numeric &n)
{
    assert(n.is_rational());
    const numeric num = n.numer();
    const numeric den = n.denom();
    assert((LONG_MIN <= num) && (num <= LONG_MAX));
    assert((LONG_MIN <= den) && (den <= LONG_MAX));
    fmpq_set_si(x, num.to_long(), den.to_long());
}

void
rat_of_ginac(rat_t rat, const GiNaC::ex &expr)
{
    LOGME;
    rat_one(rat);
    vector<ulong> exp(nvariables);
    factor_iter(expr, [&](const ex &polyfactor, int pfpower) {
        if (is_a<numeric>(polyfactor)) {
            numeric npf = ex_to<numeric>(polyfactor);
            assert(npf.is_rational());
            fmpq_t n;
            fmpq_init(n);
            fmpq_of_ginac(n, npf);
            rat_mul_fmpq(rat, n, pfpower);
            fmpq_clear(n);
        } else {
            fmpz_mpoly_t poly;
            fmpz_mpoly_init(poly, ctx);
            term_iter(polyfactor.expand(), [&](const ex &term) {
                fmpz_t coef;
                fmpz_init_set_ui(coef, 1);
                for (ulong i = 0; i < nvariables; i++) {
                    exp[i] = 0;
                }
                factor_iter(term, [&](const ex &f, int tfpower) {
                    assert(tfpower >= 0);
                    if (is_a<numeric>(f)) {
                        numeric ntf = ex_to<numeric>(f);
                        assert(ntf.is_integer());
                        fmpz_t npow;
                        fmpz_init(npow);
                        fmpz_of_ginac(npow, ntf);
                        fmpz_pow_ui(npow, npow, tfpower);
                        fmpz_mul(coef, coef, npow);
                        fmpz_clear(npow);
                    } else if (is_a<symbol>(f)) {
                        symbol sym = ex_to<symbol>(f);
                        int varidx = variable2index[sym.get_name()];
                        exp[varidx] += tfpower;
                    } else {
                        assert(!"A term has a factor which is neither a number nor a symbol");
                    }
                });
                fmpz_mpoly_push_term_fmpz_ui(poly, coef, &exp[0], ctx);
                fmpz_clear(coef);
            });
            fmpz_mpoly_sort_terms(poly, ctx);
            fmpz_mpoly_combine_like_terms(poly, ctx);
            rat_mul_fmpz_mpoly_setx(rat, poly, pfpower);
            fmpz_mpoly_clear(poly, ctx);
        }
    });
}

void
ratsum_of_ginac(ratsum_t sum, const GiNaC::ex &expr)
{
    ratsum_zero(sum);
    term_iter(expr, [&](const ex &term) {
        rat_t rat;
        rat_init(rat);
        rat_of_ginac(rat, term);
        rat_cofactorize(rat);
        ratsum_add_rat_setx(sum, rat);
        rat_clear(rat);
    });
}

/* Main
 */

int
print_ptr_arginfo(const struct printf_info *info, size_t n, int *argtypes)
{
    if (n > 0) argtypes[0] = PA_POINTER;
    return 1;
}

int
print_fmpz(FILE *f, const struct printf_info *info, const void *const *args)
{
    fmpz_fprint(f, **((const fmpz_t **)(args[0])));
    return 1;
}

int
print_poly(FILE *f, const struct printf_info *info, const void *const *args)
{
    fmpz_mpoly_fprint_pretty(f, **((const fmpz_mpoly_t **)(args[0])), variable_names, ctx);
    return 1;
}

int
print_poly_short(FILE *f, const struct printf_info *info, const void *const *args)
{
    slong nterms = fmpz_mpoly_length(**((const fmpz_mpoly_t **)(args[0])), ctx);
    fprintf(f, "P{%ldt}", nterms);
    return 1;
}

int
print_rat(FILE *f, const struct printf_info *info, const void *const *args)
{
    rat_fprint(f, **((const rat_t **)(args[0])));
    return 1;
}

int
print_rat_short(FILE *f, const struct printf_info *info, const void *const *args)
{
    rat_fprint_short(f, **((const rat_t **)(args[0])));
    return 1;
}

int
print_ratsum(FILE *f, const struct printf_info *info, const void *const *args)
{
    ratsum_fprint(f, **((const ratsum_t **)(args[0])));
    return 1;
}

GiNaC::ex
load_input(const char *filename)
{
    LOGME;
    GiNaC::parser reader;
    ex expr;
    if (strcmp(filename, "-") != 0) {
        ifstream ifs(filename);
        expr = reader(ifs);
    } else {
        expr = reader(cin);
    }
    // List variable names, init context.
    auto &name2sym = reader.get_syms();
    variable_names = (const char**)flint_calloc(name2sym.size(), sizeof(const char*));
    for (auto &&kv : name2sym) {
        int index = nvariables;
        variables.push_back(kv.first);
        variable_names[index] = (const char*)flint_calloc(variables[index].length() + 1, 1);
        strcpy((char*)variable_names[index], variables[index].c_str());
        variable2index[kv.first] = index;
        logd("Variable %d is %s", index+1, variable_names[index]);
    }
    rat_ctx_init();
    return expr;
}

void
save_output(const char *filename, rat_t rat)
{
    LOGME;
    logd("Saving %r to %s", rat, filename);
    FILE *f;
    if (strcmp(filename, "-") != 0) {
        f = fopen(filename, "wb");
    } else {
        f = stdout;
    }
    fprintf(f, "(\n");
    rat_fprint(f, rat);
    fprintf(f, "\n)\n");
    if (strcmp(filename, "-") != 0) {
        fclose(f);
    } else {
        fflush(f);
    }
}

void
usage(FILE *f)
{
    const char *p = strchr(usagetext, '\n') + 1;
    for (;;) {
        const char *l = strchr(p + 2, '{');
        if (l == NULL) break;
        const char *r = strchr(l, '}');
        if (r == NULL) break;
        const char *a = "", *b = "\033[0m";
        if (l[-2] == 'S' && l[-1] == 's') { a = "\033[1m"; }
        if (l[-2] == 'N' && l[-1] == 'm') { a = "\033[1;35m"; }
        if (l[-2] == 'F' && l[-1] == 'l') { a = "\033[33m"; }
        if (l[-2] == 'C' && l[-1] == 'm') { a = "\033[1m"; }
        if (l[-2] == 'A' && l[-1] == 'r') { a = "\033[32m"; }
        if (l[-2] == 'D' && l[-1] == 'l') { a = "\033[36m"; }
        if (l[-2] == 'Q' && l[-1] == 'l') { a = "\033[36m"; }
        fwrite(p, l - p - 2, 1, f);
        if (COLORS) fputs(a, f);
        fwrite(l + 1, r - l - 1, 1, f);
        if (COLORS) fputs(b, f);
        p = r + 1;
    }
    fputs(p, f);
}

int
main(int argc, char *argv[])
{
    register_printf_function('F', print_fmpz, print_ptr_arginfo);
    register_printf_function('P', print_poly, print_ptr_arginfo);
    register_printf_function('p', print_poly_short, print_ptr_arginfo);
    register_printf_function('R', print_rat, print_ptr_arginfo);
    register_printf_function('r', print_rat_short, print_ptr_arginfo);
    register_printf_function('S', print_ratsum, print_ptr_arginfo);
    int nthreads = 1;
    const char *inputfile = "-";
    const char *outputfile = "-";
    for (int opt; (opt = getopt(argc, (char*const*)argv, "j:hC")) != -1;) {
        switch (opt) {
        case 'j': nthreads = atoi(optarg); break;
        case 'h': usage(stdout); return 0;
        case 'C': COLORS = true; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc == 1) {
        inputfile = argv[0];
    } else if (argc == 2) {
        inputfile = argv[0];
        outputfile = argv[1];
    } else {
        fprintf(stderr, "ratnormal: bad invocation (use -h to see usage)\n");
        return 1;
    }
    flint_set_num_threads(nthreads);
    GiNaC::ex expr = load_input(inputfile);
    ratsum_t sum;
    ratsum_init(sum);
    ratsum_of_ginac(sum, expr);
    rat_t rat;
    rat_init(rat);
    ratsum_sum_setx(rat, sum);
    ratsum_clear(sum);
    rat_sort(rat);
    rat_reverse(rat);
    save_output(outputfile, rat);
    rat_clear(rat);
    rat_ctx_clear();
    return 0;
}
