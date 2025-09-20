#include <stdio.h>
#include <string.h>
#include <gmp.h>

int digit(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'z') return c - 'a' + 10;
    return -1;
}

void convert(mpz_t res, char* val, int base) {
    mpz_set_ui(res, 0);
    for (int i = 0; i < strlen(val); i++) {
        mpz_mul_ui(res, res, base);
        mpz_add_ui(res, res, digit(val[i]));
    }
}

void lagrange(mpz_t c, mpz_t y[], int x[], int k) {
    mpz_set_ui(c, 0);
    mpz_t num, den, temp;
    mpz_init(num); mpz_init(den); mpz_init(temp);
    
    for (int i = 0; i < k; i++) {
        mpz_set_ui(num, 1);
        mpz_set_ui(den, 1);
        
        for (int j = 0; j < k; j++) {
            if (i != j) {
                mpz_mul_si(num, num, -x[j]);
                mpz_mul_si(den, den, x[i] - x[j]);
            }
        }
        
        mpz_mul(temp, y[i], num);
        mpz_divexact(temp, temp, den);
        mpz_add(c, c, temp);
    }
    
    mpz_clear(num); mpz_clear(den); mpz_clear(temp);
}

int main() {
    // TEST CASE 1
    int x1[] = {1, 2, 3};
    mpz_t y1[3], c1;
    mpz_init(y1[0]); mpz_init(y1[1]); mpz_init(y1[2]); mpz_init(c1);
    
    convert(y1[0], "4", 10);
    convert(y1[1], "111", 2);
    convert(y1[2], "12", 10);
    
    lagrange(c1, y1, x1, 3);
    gmp_printf("%Zd\n", c1);
    
    // TEST CASE 2
    int x2[] = {1, 2, 3, 4, 5, 6, 7};
    mpz_t y2[7], c2;
    for (int i = 0; i < 7; i++) mpz_init(y2[i]);
    mpz_init(c2);
    
    convert(y2[0], "13444211440455345511", 6);
    convert(y2[1], "aed7015a346d635", 15);
    convert(y2[2], "6aeeb69631c227c", 15);
    convert(y2[3], "e1b5e05623d881f", 16);
    convert(y2[4], "316034514573652620673", 8);
    convert(y2[5], "2122212201122002221120200210011020220200", 3);
    convert(y2[6], "20120221122211000100210021102001201112121", 3);
    
    lagrange(c2, y2, x2, 7);
    gmp_printf("%Zd\n", c2);
    
    return 0;
}