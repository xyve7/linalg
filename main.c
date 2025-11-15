
#include <stdio.h>
#include <linalg.h>

int main() {
    Mat A = mat_from(2, 2,
        4.0, 7.0,
        2.0, 6.0
    );

    Mat Ai = mat_inverse(&A);
    Mat Aexp = mat_from(2, 2,
        0.6, -0.7,
       -0.2,  0.4
    );

    printf("Inverse of A:\n");
    mat_print(&Ai);
    printf("Expected:\n");
    mat_print(&Aexp);

    mat_free(&A);
    mat_free(&Ai);
    mat_free(&Aexp);

    Mat B = mat_from(3, 3,
        1.0, 2.0, 3.0,
        0.0, 1.0, 4.0,
        5.0, 6.0, 0.0
    );

    Mat Bi = mat_inverse(&B);
    Mat Bexp = mat_from(3, 3,
        -24.0, 18.0, 5.0,
         20.0,-15.0,-4.0,
         -5.0,  4.0, 1.0
    );

    printf("Inverse of B:\n");
    mat_print(&Bi);
    printf("Expected:\n");
    mat_print(&Bexp);

    mat_free(&B);
    mat_free(&Bi);
    mat_free(&Bexp);

    Mat C = mat_from(3, 3,
        2.0, 5.0, 3.0,
        1.0,-2.0,-1.0,
        1.0, 3.0, 4.0
    );

    Mat Ci = mat_inverse(&C);
    Mat Cexp = mat_from(3, 3,
         0.35,  0.05, -0.25,
         0.15, -0.05, -0.25,
        -0.10,  0.10,  0.20
    );

    printf("Inverse of C:\n");
    mat_print(&Ci);
    printf("Expected:\n");
    mat_print(&Cexp);

    mat_free(&C);
    mat_free(&Ci);
    mat_free(&Cexp);

    return 0;
}

