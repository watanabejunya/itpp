#include <itpp/itbase.h>
#include <itpp/comm/llr.h>
#include "../../Library/Common/helper.hpp"
#include "../../Library/TurboCode/RecursiveConvolutionalCode.hpp"

using std::max;
using std::string;
using itpp::bin;
using itpp::bvec;
using itpp::ivec;
using itpp::vec;
using itpp::max;
using itpp::pow2i;
using itpp::sum;
using itpp::elem_mult;
using itpp::dec2bin;
using itpp::bin2dec;
using itpp::reverse_int;
using itpp::LLR_calc_unit;
using helper::INFTY;


RecursiveConvolutionalCode::RecursiveConvolutionalCode()
{
    this->LLRCalculator = LLR_calc_unit();
    this->set_termination_type(true);
}

RecursiveConvolutionalCode::~RecursiveConvolutionalCode()
{
}

void RecursiveConvolutionalCode::set_generator_polynomials(const ivec& polynomial, const int constraintLength)
{
    it_assert_debug(polynomial.size() == 2, "RecursiveConvolutionalCode::set_generator_polynomials(): Polynomial size must be 2.");
    it_assert_debug(polynomial[0] < pow2i(constraintLength), "RecursiveConvolutionalCode::set_generator_polynomials(): Too small constraint length.");
    it_assert_debug(polynomial[1] < pow2i(constraintLength), "RecursiveConvolutionalCode::set_generator_polynomials(): Too small constraint length.");

    this->backwardPolynomial = dec2bin(constraintLength, polynomial[0]);
    this->forwardPolynomial = dec2bin(constraintLength, polynomial[1]);
    this->constraintLength = constraintLength;
    this->memorySize = constraintLength - 1;
    this->numberOfStates = pow2i(this->memorySize);
    this->outputBlockLength = polynomial.size();
    this->rate = 1.0 / polynomial.size();

    this->stateTable.set_size(this->numberOfStates, 2);
    this->outputTable.set_size(this->numberOfStates, 2 * (this->outputBlockLength - 1));

    // Calculate next state and output and set on the table.
    for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
        for (int input = 0; input <= 1; input++) {
            this->stateTable(currentState, input) = calculate_next_state(currentState, input);
            this->outputTable(currentState, input) = calculate_output(currentState, input);
        }
    }
}

void RecursiveConvolutionalCode::set_method_type(const string& method)
{
    this->methodType = method;
}

void RecursiveConvolutionalCode::set_termination_type(const bool terminated)
{
    this->isTerminated = terminated;
}

int RecursiveConvolutionalCode::calculate_output(const int state, const int input)
{
    bvec memories;

    it_assert_debug(0 <= state && state < this->numberOfStates, "RecursiveConvolutionalCode::calculate_output(): Invalid state given.");
    it_assert_debug(input == 0 || input == 1, "RecursiveConvolutionalCode::calculate_output(): Invalid input given.");

    // Set state and input on memories.
    memories = dec2bin(this->constraintLength, state);
    memories[0] = bin(input);

    // Calculate the feedback bit by XOR-addition.
    memories[0] = sum(elem_mult(memories, this->backwardPolynomial));

    // Calculate the feedforward bit by XOR-addition.
    return sum(elem_mult(memories, this->forwardPolynomial));
}

int RecursiveConvolutionalCode::calculate_next_state(const int state, const int input)
{
    bvec memories;

    it_assert_debug(0 <= state && state < this->numberOfStates, "RecursiveConvolutionalCode::calculate_next_state(): Invalid state given.");
    it_assert_debug(input == 0 || input == 1, "RecursiveConvolutionalCode::calculate_next_state(): Invalid input given.");

    // Set state and input on memories.
    memories = dec2bin(this->constraintLength, state);
    memories[0] = bin(input);

    // Calculate the feedback bit by XOR-addition.
    memories[0] = sum(elem_mult(memories, this->backwardPolynomial));

    // Convert the binary vector to an integer that presents next state.
    return bin2dec(memories.left(this->memorySize));
}

double RecursiveConvolutionalCode::get_rate() const
{
    return this->rate;
}

bvec RecursiveConvolutionalCode::encode(const bvec& input)
{
    bvec output;
    this->encode(input, output);
    return output;
}

void RecursiveConvolutionalCode::encode(const bvec& input, bvec& output)
{
    if (this->isTerminated) {
        this->encode_tail(input, output);
    } else {
        this->encode_trunc(input, output);
    }
}

void RecursiveConvolutionalCode::encode_trunc(const bvec& input, bvec& output)
{
    const int informationLength = input.size();
    int state = 0;

    output.set_size(informationLength);

    for (int i = 0; i < informationLength; i++) {
        output[i] = this->outputTable(state, int(input[i]));
        state = this->stateTable(state, int(input[i]));
    }
}

void RecursiveConvolutionalCode::encode_tail(const bvec& input, bvec& output)
{
    const int informationLength = input.size();
    int state = 0, targetState;

    output.set_size(informationLength + 2 * this->memorySize);

    for (int i = 0; i < informationLength; i++) {
        output[i] = this->outputTable(state, int(input[i]));
        state = this->stateTable(state, int(input[i]));
    }

    // Add K-1 zeros as tail bits.
    for (int i = informationLength; i < informationLength + this->memorySize; i++) {
        targetState = (state << 1) & (pow2i(this->memorySize) - 1);
        output[i] = this->stateTable(state, 0) == targetState ? bin(0) : bin(1);
        output[i + 1] = this->outputTable(state, int(output[i]));
        state = targetState;
    }
}

void RecursiveConvolutionalCode::decode(const vec& systematicBits, const vec& parityBits, const vec& extrinsicLLRs, vec& LLRs)
{
    if (this->methodType == "TABLE") {
        this->decode_table(systematicBits, parityBits, extrinsicLLRs, LLRs);
    } else if (this->methodType == "MAP") {
        this->decode_map(systematicBits, parityBits, extrinsicLLRs, LLRs);
    } else if (this->methodType == "LOGMAP") {
        this->decode_log_map(systematicBits, parityBits, extrinsicLLRs, LLRs);
    } else if (this->methodType == "LOGMAX") {
        this->decode_log_max(systematicBits, parityBits, extrinsicLLRs, LLRs);
    }
}

// void RecursiveConvolutionalCode::decode_map(const vec& rec_systematic, const mat &rec_parity, const vec& extrinsicLLRs, vec& LLRs)
// {
//     int codeLength;
//     double gamma_k_e, nom, den, temp0, temp1, exp_temp0, exp_temp1;
//     int j, s0, s1, node, node - 1, s, s_prim0, s_prim1, inputLength = rec_systematic.length();
//     ivec p0, p1;
//
//     mat alpha(this->numberOfStates, inputLength + 1);
//     mat beta(this->numberOfStates, inputLength + 1);
//     mat gamma(2 * this->numberOfStates, inputLength + 1);
//     vec denom(inputLength + 1);
//     denom.clear();
//
//     LLRs.set_size(inputLength, false);
//
//     //Calculate gamma
//     for (int node = 1; node <= inputLength; node++) {
//         for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
//             exp_temp0 = 0.0;
//             exp_temp1 = 0.0;
//             for (j = 0; j < (this->outputBlockLength - 1); j++) {
//                 exp_temp0 += 0.5 * rec_parity(node - 1, j) * double(1 - 2 * this->outputTable(currentState, 2 * j + 0)); /* Is this OK? */
//                 exp_temp1 += 0.5 * rec_parity(node - 1, j) * double(1 - 2 * this->outputTable(currentState, 2 * j + 1));
//             }
//
//             gamma(2*currentState + 0, node) = trunc_exp(0.5 * (extrinsicLLRs(node - 1) + Lc * rec_systematic(node - 1)) + exp_temp0);
//             gamma(2*currentState + 1, node) = trunc_exp(-0.5 * (extrinsicLLRs(node - 1) + Lc * rec_systematic(node - 1)) + exp_temp1);
//         }
//     }
//
//     //Initiate alpha
//     alpha.set_col(0, zeros(this->numberOfStates));
//     alpha(0, 0) = 1.0;
//
//     //Calculate alpha and denom going forward through the trellis
//     for (node = 1; node <= inputLength; node++) {
//         for (s = 0; s < this->numberOfStates; s++) {
//             s_prim0 = this->reverseStateTable(s, 0);
//             s_prim1 = this->reverseStateTable(s, 1);
//             temp0 = alpha(s_prim0, node - 1) * gamma(2 * s_prim0 + 0, node);
//             temp1 = alpha(s_prim1, node - 1) * gamma(2 * s_prim1 + 1, node);
//             alpha(s, node) = temp0 + temp1;
//             denom(node)  += temp0 + temp1;
//         }
//         alpha.set_col(node, alpha.get_col(node) / denom(node));
//     }
//
//     //Initiate beta
//     if (this->isTerminated) {
//         beta.set_col(inputLength, zeros(this->numberOfStates));
//         beta(0, inputLength) = 1.0;
//     } else {
//         beta.set_col(inputLength, alpha.get_col(inputLength));
//     }
//
//     //Calculate beta going backward in the trellis
//     for (node = inputLength; node >= 2; node--) {
//         for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//             s0 = this->stateTable(currentState, 0);
//             s1 = this->stateTable(currentState, 1);
//             beta(currentState, node - 1) = (beta(s0, node) * gamma(2 * currentState + 0, node)) + (beta(s1, node) * gamma(2 * currentState + 1, node));
//         }
//         beta.set_col(node - 1, beta.get_col(node - 1) / denom(node));
//     }
//
//     //Calculate extrinsic output for each bit
//     for (node = 1; node <= inputLength; node++) {
//         node - 1 = node - 1;
//         nom = 0;
//         den = 0;
//         for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//             s0 = this->stateTable(currentState, 0);
//             s1 = this->stateTable(currentState, 1);
//             exp_temp0 = 0.0;
//             exp_temp1 = 0.0;
//             exp_temp0 += 0.5 * Lc * rec_parity(node - 1, j) * double(1 - 2 * this->outputTable(currentState, 2 * j + 0));
//             exp_temp1 += 0.5 * Lc * rec_parity(node - 1, j) * double(1 - 2 * this->outputTable(currentState, 2 * j + 1));
//             // gamma_k_e = std::exp( exp_temp0 );
//             gamma_k_e = trunc_exp(exp_temp0);
//             nom += alpha(currentState, node - 1) * gamma_k_e * beta(s0, node);
//
//             // gamma_k_e = std::exp( exp_temp1 );
//             gamma_k_e = trunc_exp(exp_temp1);
//             den += alpha(currentState, node - 1) * gamma_k_e * beta(s1, node);
//         }
//         //      LLRs(node - 1) = std::log(nom/den);
//         LLRs(node - 1) = trunc_log(nom / den);
//     }
//
// }

void RecursiveConvolutionalCode::decode_log_max(const vec& rec_systematic, const vec& rec_parity,
                                       const vec& extrinsicLLRs, vec& LLRs)
{
    const int inputLength = rec_systematic.size(), LLRLength = extrinsicLLRs.size();
    double nom, den, exp_temp, ex, rp;

    mat alpha(this->numberOfStates, inputLength + 1);
    mat beta(this->numberOfStates, inputLength + 1);
    mat gamma(2 * this->numberOfStates, inputLength + 1);
    LLRs.set_size(LLRLength, false);

    it_assert_debug(inputLength <= LLRLength, "RecursiveConvolutionalCode::decode_log_max(): Too short LLR length.");

    // Initiate alpha.
    alpha(0, 0) = 0;
    alpha.set_submatrix(1, this->numberOfStates - 1, 0, 0, -INFTY);

    // Calculate alpha and gamma going forward through the trellis.
    for (int node = 1; node <= inputLength; node++) {
        ex = (extrinsicLLRs[node - 1] + rec_systematic[node - 1]) / 2.0;
        rp = rec_parity[node - 1] / 2.0;

        for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
            double metric = -INFTY;
            for (int input = 0; input <= 1; input++) {
                previousState = this->reverseStateTable(currentState, input);
                exp_temp = (this->outputTable(previousState, input) == 1) ? -rp : rp;
                gamma[2 * previousState + input, node] = (input) ? exp_temp + ex : exp_temp - ex;

                metric = max(metric, alpha(previousState, node - 1) + gamma(2 * previousState + input, node));
            }
            alpha(currentState, node) = metric;
        }

        // Normalize metric.
        alpha.set_column(node) = alpha.get_column(node) - alpha(0, node);
    }

    // Initiate beta.
    if (this->isTerminated) {
        beta(0, inputLength) = 0;
        beta.set_submatrix(1, this->numberOfStates - 1, inputLength, inputLength, -INFTY);
    } else {
        beta.set_col(inputLength, alpha.get_col(inputLength));
    }

    // Calculate beta going backward in the trellis.
    for (int node = inputLength; node >= 1; node--) {
        for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
            double metric = -INFTY;
            for (int input = 0; input <= 1; input++) {
                metric = max(metric, alpha(currentState, node - 1) + gamma(2 * currentState + input, node));
                if (metric < beta(this->stateTable(currentState, input), node) + gamma(2 * currentState + input, node)) {
                    metric = beta(this->stateTable(currentState, input), node) + gamma(2 * currentState + input, node)
                }
            }
            beta(currentState, node - 1) = metric;
        }

        // Normalize metric.
        beta.set_column(node) = beta.get_column(node) - beta(0, node);
    }

    // Calculate extrinsic output for each bit.
    for (int node = 1; node <= inputLength; node++) {
        double metricFor0 = -INFTY;
        double metricFor1 = -INFTY;
        for (int currentState = 0; currentState < this->numberOfStates; currentState++) {
            for (int input = 0; input <= 1; input++) {
                int nextState = this->stateTable(currentState, input);
                int output = this->outputTable(currentState , input);
                double gamma = (output == 0) ? rec_parity[node - 1] / 2.0 : - rec_parity[node - 1] / 2.0;

                if (output == 0) {
                    metricFor0 = max(metricFor0, alpha(currentState, node - 1) + gamma + beta(nextState, node));
                } else if (output == 1) {
                    metricFor1 = max(metricFor1, alpha(currentState, node - 1) + gamma + beta(nextState, node));
                }
            }
        }
        LLRs[node - 1] = metricFor0 - metricFor1;
    }
}

// void RecursiveConvolutionalCode::decode_log_map(const vec& rec_systematic, const vec& rec_parity,
//                                        const vec& extrinsicLLRs, vec& LLRs)
// {
//     const int inputLength = rec_systematic.size();
//
//     double nom, den, exp_temp0, exp_temp1, rp;
//     int node, node - 1, l, s, currentState, s_prim0, s_prim1;
//     int LLRLength = extrinsicLLRs.length();
//     ivec p0, p1;
//     double ex, norm;
//     mat alpha(this->numberOfStates, inputLength + 1);
//     mat beta(this->numberOfStates, inputLength + 1);
//     mat gamma(2 * this->numberOfStates, inputLength + 1);
//
//     LLRs.set_size(LLRLength, false);
//     //denom.set_size(inputLength+1,false); for (node=0; node<=inputLength; node++) { denom(node) = -INFTY; }
//
//     //Initiate alpha
//     alpha(0, 0) = 0;
//     alpha.set_submatrix(1, this->numberOfStates - 1, 0, 0, -INFTY);
//
//     //Calculate alpha and gamma going forward through the trellis
//     for (int node = 1; node <= inputLength; node++) {
//         if (node - 1 < LLRLength) {
//             ex = 0.5 * (extrinsicLLRs(node - 1) + rec_systematic(node - 1));
//         } else {
//             ex = 0.5 * rec_systematic(node - 1);
//         }
//         rp = 0.5 * rec_parity(node - 1);
//         for (s = 0; s < this->numberOfStates; s++) {
//         s_prim0 = this->reverseStateTable(s, 0);
//         s_prim1 = this->reverseStateTable(s, 1);
//         if (this->outputTable(s_prim0 , 0)) { exp_temp0 = -rp; }
//         else { exp_temp0 = rp; }
//         if (this->outputTable(s_prim1 , 1)) { exp_temp1 = -rp; }
//         else { exp_temp1 = rp; }
//         gamma(2*s_prim0  , node) =   ex + exp_temp0;
//         gamma(2*s_prim1 + 1, node) =  -ex + exp_temp1;
//         alpha(s, node) = log_add(alpha(s_prim0, node - 1) + gamma(2 * s_prim0  , node), alpha(s_prim1, node - 1) + gamma(2 * s_prim1 + 1, node));
//         //denom(node)   = log_add( alpha(s,node), denom(node) );
//     }
//         norm = alpha(0, node); //norm = denom(node);
//         for (l = 0; l < this->numberOfStates; l++) { alpha(l, node) -= norm; }
//     }
//
//     //Initiate beta
//     if (this->isTerminated) {
//         for (s = 1; s < this->numberOfStates; s++) {
//             beta(s, inputLength) = -INFTY;
//         }
//         beta(0, inputLength) = 0.0;
//     } else {
//         beta.set_col(inputLength, alpha.get_col(inputLength));
//     }
//
//     //Calculate beta going backward in the trellis
//     for (node = inputLength; node >= 1; node--) {
//         node - 1 = node - 1;
//         for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//             beta(currentState, node - 1) = log_add(beta(this->stateTable(currentState, 0), node) + gamma(2 * currentState, node),
//                                  beta(this->stateTable(currentState, 1), node) + gamma(2 * currentState + 1, node));
//         }
//         norm = beta(0, node); //norm = denom(node);
//         for (l = 0; l < this->numberOfStates; l++) { beta(l, node) -= norm; }
//     }
//
//     //Calculate extrinsic output for each bit
//     for (node = 1; node <= inputLength; node++) {
//         node - 1 = node - 1;
//         if (node - 1 < LLRLength) {
//             nom = -INFTY;
//             den = -INFTY;
//             rp = 0.5 * rec_parity(node - 1);
//             for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//                 if (this->outputTable(currentState , 0)) {
//                     exp_temp0 = -rp;
//                 } else {
//                     exp_temp0 = rp;
//                 }
//                 if (this->outputTable(currentState , 1)) {
//                     exp_temp1 = -rp;
//                 } else {
//                     exp_temp1 = rp;
//                 }
//                 nom = log_add(nom, alpha(currentState, node - 1) + exp_temp0 + beta(this->stateTable(currentState, 0), node));
//                 den = log_add(den, alpha(currentState, node - 1) + exp_temp1 + beta(this->stateTable(currentState, 1), node));
//             }
//             LLRs(node - 1) = nom - den;
//         }
//     }
//
// }
//
// void RecursiveConvolutionalCode::table_decode(const vec&rec_systematic,
//                                        const vec&rec_parity,
//                                        const vec&extrinsicLLRs,
//                                        vec& LLRs)
// {
//   int nom, den, exp_temp0, exp_temp1, rp;
//   int node, node - 1, l, s, currentState, s_prim0, s_prim1, inputLength = rec_systematic.length();
//   int LLRLength = extrinsicLLRs.length();
//   ivec p0, p1;
//   int ex, norm;
//
//     QLLRvec rec_systematic_q  = this->LLRCalculator.to_qllr(rec_systematic);
//     QLLRvec rec_parity_q = this->LLRCalculator.to_qllr(rec_parity);
//     QLLRvec extrinsic_input_q = this->LLRCalculator.to_qllr(extrinsicLLRs);
//     QLLRvec extrinsic_output_q(length(LLRs));
//     log_decode_n2(rec_systematic_q, rec_parity_q, extrinsic_input_q,
//                   extrinsic_output_q);
//     LLRs = this->LLRCalculator.to_double(extrinsic_output_q);
//
//   QLLRmat alpha_q(this->numberOfStates, inputLength + 1);
//   QLLRmat beta_q(this->numberOfStates, inputLength + 1);
//   QLLRmat gamma_q(2*this->numberOfStates, inputLength + 1);
//   LLRs.set_size(LLRLength, false);
//   //denom.set_size(inputLength+1,false); for (node=0; node<=inputLength; node++) { denom(node) = -INFTY; }
//
//   //Initiate alpha
//   for (s = 1; s < this->numberOfStates; s++) { alpha_q(s, 0) = -QLLR_MAX; }
//   alpha_q(0, 0) = 0;
//
//   //Calculate alpha and gamma going forward through the trellis
//   for (node = 1; node <= inputLength; node++) {
//     node - 1 = node - 1;
//     if (node - 1 < LLRLength) {
//       ex = (extrinsicLLRs(node - 1) + rec_systematic(node - 1)) / 2;
//     }
//     else {
//       ex =  rec_systematic(node - 1) / 2;
//     }
//     rp =  rec_parity(node - 1) / 2;
//     for (s = 0; s < this->numberOfStates; s++) {
//       s_prim0 = this->reverseStateTable(s, 0);
//       s_prim1 = this->reverseStateTable(s, 1);
//       if (this->outputTable(s_prim0 , 0)) { exp_temp0 = -rp; }
//       else { exp_temp0 = rp; }
//       if (this->outputTable(s_prim1 , 1)) { exp_temp1 = -rp; }
//       else { exp_temp1 = rp; }
//       gamma_q(2*s_prim0  , node) =   ex + exp_temp0;
//       gamma_q(2*s_prim1 + 1, node) =  -ex + exp_temp1;
//       alpha_q(s, node) = this->LLRCalculator.jaclog(alpha_q(s_prim0, node - 1) + gamma_q(2 * s_prim0  , node),
//                                      alpha_q(s_prim1, node - 1) + gamma_q(2 * s_prim1 + 1, node));
//       //denom(node)   = log_add( alpha(s,node), denom(node) );
//     }
//     norm = alpha_q(0, node); //norm = denom(node);
//     for (l = 0; l < this->numberOfStates; l++) { alpha_q(l, node) -= norm; }
//   }
//
//   //Initiate beta
//   if (this->isTerminated) {
//     for (s = 1; s < this->numberOfStates; s++) { beta_q(s, inputLength) = -QLLR_MAX; }
//     beta_q(0, inputLength) = 0;
//   }
//   else {
//     beta_q.set_col(inputLength, alpha_q.get_col(inputLength));
//   }
//
//   //Calculate beta going backward in the trellis
//   for (node = inputLength; node >= 1; node--) {
//     node - 1 = node - 1;
//     for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//       beta_q(currentState, node - 1) = this->LLRCalculator.jaclog(beta_q(this->stateTable(currentState, 0), node) + gamma_q(2 * currentState, node),
//                                           beta_q(this->stateTable(currentState, 1), node) + gamma_q(2 * currentState + 1, node));
//     }
//     norm = beta_q(0, node); //norm = denom(node);
//     for (l = 0; l < this->numberOfStates; l++) { beta_q(l, node) -= norm; }
//   }
//
//   //Calculate extrinsic output for each bit
//   for (node = 1; node <= inputLength; node++) {
//     node - 1 = node - 1;
//     if (node - 1 < LLRLength) {
//       nom = -QLLR_MAX;
//       den = -QLLR_MAX;
//       rp =  rec_parity(node - 1) / 2;
//       for (currentState = 0; currentState < this->numberOfStates; currentState++) {
//         if (this->outputTable(currentState , 0)) { exp_temp0 = -rp; }
//         else { exp_temp0 = rp; }
//         if (this->outputTable(currentState , 1)) { exp_temp1 = -rp; }
//         else { exp_temp1 = rp; }
//         nom = this->LLRCalculator.jaclog(nom, alpha_q(currentState, node - 1) + exp_temp0 + beta_q(this->stateTable(currentState, 0), node));
//         den = this->LLRCalculator.jaclog(den, alpha_q(currentState, node - 1) + exp_temp1 + beta_q(this->stateTable(currentState, 1), node));
//       }
//       LLRs(node - 1) = nom - den;
//     }
//   }
// }
