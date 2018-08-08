#include <cstring>
#include <iostream>
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "parasail/matrix_lookup.h"

int main(int argc, char **argv) {
    // align two sequences s1, and s2
    const char *s1 = "asdfasdfffffffasdf";
    const char *s2 = "asdfasdfkffbffasdf";
    // sequence lengths needed for parasail function call
    int s1Len = (int)strlen(s1);
    int s2Len = (int)strlen(s2);

    // initialize parasail result
    parasail_result_t *result = NULL;
    // initialize parasail scoring matrix
    const parasail_matrix_t *matrix = NULL;

    /* note the address-of operator '&' */
    // parasail_blosum62, a score matrix pre-defined in "parasail/matrices/blosum62.h"
    // perform smith-waterman local alignment, no trace
    result = parasail_sw(s1, s1Len, s2, s2Len, 11, 1, &parasail_blosum62);

    // free result
    parasail_result_free(result);

    // use another score matrix "pam100"
    matrix = parasail_matrix_lookup("pam100");
    // perform sw local alignment with traces
    result = parasail_sw_trace(s1, s1Len, s2, s2Len, 11, 1, matrix);

    // generic trackback function to visualize alignments
    // s1 name="Query:", s2 name="Target:", '|' for match and '*' for mismatches
    // width=60
    parasail_traceback_generic(s1, s1Len, s2, s2Len, "Query:", "Target:", matrix, result, '|', '*', '*', 60, 7, 1);
    parasail_cigar_* cigar = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, matrix);

    // print out cigar length, query s1 start pos, target s2 start pos
    std::cout << std::endl << "print cigar " << cigar->len << ", " << cigar->beg_query << ", " << cigar->beg_ref << std::endl;
    for (auto i=0; i<cigar->len; ++i) {
        std::cout << parasail_cigar_decode_len(cigar->seq[i])  // cigar length
                  << parasail_cigar_decode_op(cigar->seq[i]) << " ";  // cigar op
    }
    std::cout << std::endl;

    // free result
    parasail_result_free(result);
    // free cigar
    parasail_cigar_free(cigar);

// Expected Output:
//                 Target:       1 asdfasdfkffbffasdf      18
//                 ||||||||*||*||||||
//                 Query:        1 asdfasdfffffffasdf      18
//
//                 Length: 18
//                 Identity:        16/18 (88.9%)
//                 Similarity:      18/18 (100.0%)
//                 Gaps:             0/18 ( 0.0%)
//                 Score: 82
//
//                 print cigar 5, 0, 0
//                 8= 1X 2= 1X 6=
}