

#include <stdio.h>
#include <stdlib.h>


#include <R.h>
#include <Rinternals.h>


SEXP hermite(SEXP n)
{
    SEXP value = PROTECT(allocVector(REALSXP, 3));
    UNPROTECT(1);
    return value;
}


//double * hermite(int n)
//{
//    if (n < 0) {
//        printf("'n' must be non-negative\n");
//        double value[] = {1.0};
//        return &value;
//    }
//    double value[n + 1];
//    value[0] = 1.0;
//    for (int i = 0; i < n; i++) {
//        double tmp[i];
//        for (int j = 0; j < i; j++)
//            tmp[j] = value[j] - value[j + 2] * (2.0 + (double) j);
//        value[0] = -value[1];
//        value[i + 1] = 1.0;
//        for (int j = 0; j < i; j++)
//            value[j + 1] = tmp[j];
//    }
//    for (int k = 0; k < n + 1; k++)
//        printf("%f\n", value[k]);
//    return &value;
//}


//int main()
//{
//    hermite(10);
//    return 0;
//}
