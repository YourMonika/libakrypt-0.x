/* ----------------------------------------------------------------------------------------------- */
/*      Реализация блочного алгоритма шифрования TwoFish (128-bit), author: Bruce Schneier         */
/* ----------------------------------------------------------------------------------------------- */

#include <math.h>
#include <libakrypt.h>
#include <assert.h>


#define RHO 0x01010101L
/*----------------------------------------------------------------
 *  FUNCTIONAL STRUCT
 -----------------------------------------------------------------*/
typedef struct {
    UINT Kkey[40];
    UINT Skey[2];
} ak_key_shedule;

/*----------------------------------------------------------------
 *  CIRCULAR SHIFTS
 -----------------------------------------------------------------*/
#define ROL(x,n) (((x) << ((n) & 0x1F)) | ((x) >> (32-((n) & 0x1F))))
#define ROR(x,n) (((x) >> ((n) & 0x1F)) | ((x) << (32-((n) & 0x1F))))
//Right circular shift for 4 least significant bits
ak_uint8 ROR4(ak_uint8 x){
    return (((x << 3) & 0xF) | ( (x & 0xF) >> 1));
}
/*----------------------------------------------------------------
 *  CONST MATRIX 
 -----------------------------------------------------------------*/
const ak_uint8 MDS[4][4] = {
        {0x01, 0xEF, 0x5B, 0x5B},
        {0x5B, 0xEF, 0xEF, 0x01},
        {0xEF, 0x5B, 0x01, 0xEF},
        {0xEF, 0x01, 0xEF, 0x5B}
};
const ak_uint8 qt[2][4][16] = {
        //permutation table for q0 operation
        {
                { 0x8, 0x1, 0x7, 0xD, 0x6, 0xF, 0x3, 0x2, 0x0, 0xB, 0x5, 0x9, 0xE, 0xC, 0xA, 0x4 },
                { 0xE, 0xC, 0xB, 0x8, 0x1, 0x2, 0x3, 0x5, 0xF, 0x4, 0xA, 0x6, 0x7, 0x0, 0x9, 0xD },
                { 0xB, 0xA, 0x5, 0xE, 0x6, 0xD, 0x9, 0x0, 0xC, 0x8, 0xF, 0x3, 0x2, 0x4, 0x7, 0x1 },
                { 0xD, 0x7, 0xF, 0x4, 0x1, 0x2, 0x6, 0xE, 0x9, 0xB, 0x3, 0x0, 0x8, 0x5, 0xC, 0xA }
        },

        //permutation table for q1 operation
        {
                { 0x2, 0x8, 0xB, 0xD, 0xF, 0x7, 0x6, 0xE, 0x3, 0x1, 0x9, 0x4, 0x0, 0xA, 0xC, 0x5 },
                { 0x1, 0xE, 0x2, 0xB, 0x4, 0xC, 0x3, 0x7, 0x6, 0xD, 0xA, 0x5, 0xF, 0x9, 0x0, 0x8 },
                { 0x4, 0xC, 0x7, 0x5, 0x1, 0x6, 0x9, 0xA, 0x0, 0xE, 0xD, 0x8, 0x2, 0xB, 0x3, 0xF },
                { 0xB, 0x9, 0x5, 0x1, 0xC, 0x3, 0xD, 0xE, 0x6, 0x4, 0x7, 0xF, 0x2, 0x0, 0x8, 0xA }
        },
};
const ak_uint8 RSMatrix[4][8] = {
        {0x01, 0xA4, 0x55, 0x87, 0x5A, 0x58, 0xDB, 0x9E},
        {0xA4, 0x56, 0x82, 0xF3, 0x1E, 0xC6, 0x68, 0xE5},
        {0x02, 0xA1, 0xFC, 0xC1, 0x47, 0xAE, 0x3D, 0x19},
        {0xA4, 0x55, 0x87, 0x5A, 0x58, 0xDB, 0x9E, 0x03},
};
/*----------------------------------------------------------------
 *  FUNCTIONS FOR MUTATION TYPES ak_uint8, ak_uint32
 -----------------------------------------------------------------*/
ak_uint8 ak_b0(ak_uint32 x){
    return x & 0xFF;
}
ak_uint8 ak_b1(ak_uint32 x){
    return  (x >> 8) & 0xFF;
}
ak_uint8  ak_b2(ak_uint32 x){
    return  (x >> 16) & 0xFF;
}
ak_uint8  ak_b3(ak_uint32 x){
    return  (x >> 24) & 0xFF;
}
ak_uint32 ak_uint8_to_uint32(ak_uint8 x0, ak_uint8 x1, ak_uint8 x2, ak_uint8 x3){
    return ((x0 << 24) ^ (x1 << 16) ^ (x2 << 8) ^ x3);
}
ak_uint32 ak_auint8_to_uint32(const ak_uint8 *x){
    return ((x[0] << 24) ^ (x[1] << 16) ^ (x[2] << 8) ^ x[3]);
}

/*----------------------------------------------------------------
 *  FUNCTIONS FOR POLYNOMIAL MULTIPLYING, MODing, MulMODing
 -----------------------------------------------------------------*/
ak_uint32 polyMul(ak_uint32 A, ak_uint32 B){
    ak_uint32 t=0;
    while (A)
    {
        if (A&1) t^=B;
        B <<= 1;
        A >>= 1;
    }
    return t;
}

ak_uint32 gfMOD(ak_uint32 A, ak_uint32 modulus){
    ak_uint32 tt;

    modulus <<= 7;
    for (int i = 0; i < 8; i++)
    {
        tt = A ^ modulus;
        if (tt < A) A = tt;
        modulus >>= 1;
    }
    return A;
}

#define     gfMul(a, b, modulus)    gfMOD(polyMul(a, b), modulus)
/*----------------------------------------------------------------
 *             TwoFish - IN-ROUND FUNCTIONS
 -----------------------------------------------------------------*/
ak_uint8 q(ak_uint8 x, int op){
    assert( (op == 0) || (op == 1));
    //splitting ak_uint8 into two nibbles
    ak_uint8 a0 = x / 16;
    ak_uint8 b0 = x % 16;

    ak_uint8 a1 = a0 ^ b0;
    ak_uint8 b1 = a0 ^ ROR4(b0) ^ ((8 * a0) % 16); //ROR4

    ak_uint8 a2 = qt[op][0][a1];
    ak_uint8 b2 = qt[op][1][b1];

    ak_uint8 a3 = a2 ^ b2;
    ak_uint8 b3 = a2 ^ ROR4(b2) ^ ((8 * a2) % 16);  //ROR4
    ak_uint8 a4 = qt[op][2][a3];
    ak_uint8 b4 = qt[op][3][b3];

    return (16*b4 + a4);
}

ak_uint32 g(ak_uint32 X, ak_uint32 S[2]){
    ak_uint8 x0, x1, x2, x3;
    x0 = ak_b0(X);
    x1 = ak_b1(X);
    x2 = ak_b2(X);
    x3 = ak_b3(X);

    x0 = q(q(q(x0,0)^ak_b0(S[0]), 0)^ak_b0(S[1]), 1);
    x1 = q(q(q(x1,1)^ak_b1(S[0]), 0)^ak_b1(S[1]), 0);
    x2 = q(q(q(x2,0)^ak_b2(S[0]), 1)^ak_b2(S[1]), 1);
    x3 = q(q(q(x3,1)^ak_b3(S[0]), 1)^ak_b3(S[1]), 0);

    ak_uint8 x[4] = {x1, x2, x3, x0};

    ak_uint8 t;
    ak_uint8 result[4];

    for (int j = 0; j < 4; j++)
    {
        t = 0;
        for (int k = 0; k < 4; k++)
        {
            t ^= gfMul(MDS[j][k], x[k], 0x169);
        }
        result[3-j] = t;
    }
    return ak_auint8_to_uint32(result);
}


/*----------------------------------------------------------------
 *             TwoFish - ROUND FUNCTIONS
 -----------------------------------------------------------------*/
void ROUNDENCRYPT(ak_uint32 R[4], int round, ak_key_shedule keys){
    ak_uint32 T0, T1;
    T0 = g(R[0], keys.Skey);
    T1 = g(ROL(R[1], 8), keys.Skey);

    R[2] = R[0];
    R[3] = R[1];
    //adding keys from keyShedule
    R[0] = ROR(R[2] ^ (T1 + T0 + keys.Kkey[2*round + 8]), 1);
    R[1] = ROL(R[3], 1) ^ (2*T1 + T0 + keys.Kkey[2*round + 9]);

}

void ROUNDDECRYPT(ak_uint32 R[4], int round, ak_key_shedule keys){
    ak_uint32 T0, T1;
    T0 = g(R[0], keys.Skey);
    T1 = g(ROL(R[1], 8), keys.Skey);

    R[2] = R[0];
    R[3] = R[1];
    //adding keys from keyShedule
    R[0] = ROL(R[2], 1) ^ (T1 + T0 + keys.Kkey[2*round + 8]);
    R[1] = ROR(R[3] ^ (2*T1 + T0 + keys.Kkey[2*round + 9]), 1);
}

void ak_twofish_encrypt(const ak_uint32 plain[4], ak_key_shedule keys, ak_uint32 R[4]){
    R[0] = plain[0]^keys.Kkey[0];
    R[1] = plain[1]^keys.Kkey[1];
    R[2] = plain[2]^keys.Kkey[2];
    R[3] = plain[3]^keys.Kkey[3];

    for(int i=0; i<16; i++){
        ROUNDENCRYPT(R, i, keys);
    }

    R[0] = R[0]^keys.Kkey[4];
    R[1] = R[1]^keys.Kkey[5];
    R[2] = R[2]^keys.Kkey[6];
    R[3] = R[3]^keys.Kkey[7];
}

void ak_twofish_decrypt(const ak_uint32 chipher[4], ak_key_shedule keys, ak_uint32 R[4]){
    R[0] = chipher[0]^keys.Kkey[4];
    R[1] = chipher[1]^keys.Kkey[5];
    R[2] = chipher[2]^keys.Kkey[6];
    R[3] = chipher[3]^keys.Kkey[7];

    for(int i=15; i>=0; i--){
        ROUNDDECRYPT(R, i, keys);
    }

    R[0] = R[0]^keys.Kkey[0];
    R[1] = R[1]^keys.Kkey[1];
    R[2] = R[2]^keys.Kkey[2];
    R[3] = R[3]^keys.Kkey[3];
}



/*----------------------------------------------------------------
 *             TwoFish - MATRIX-MUL FUNCTIONS
 -----------------------------------------------------------------*/
ak_uint32 RSM(ak_uint8 sd[8])
{
    int j, k;
    ak_uint8 t;
    ak_uint8 result[4];

    for (j = 0; j < 4; j++)
    {
        t = 0;
        for (k = 0; k < 8; k++)
        {
            t ^= gfMul(RSMatrix[j][k], sd[k], 0x14D);
        }
        result[3-j] = t;
    }
    return ak_auint8_to_uint32(result);
}

/*----------------------------------------------------------------
 *             TwoFish - KEY-GEN FUNCTIONS
 -----------------------------------------------------------------*/
void pushSkeys(const ak_uint8 globalkey[16], ak_key_shedule key){
    ak_uint8 m07[8] = {globalkey[0], globalkey[1], globalkey[2], globalkey[3],
                   globalkey[4], globalkey[5], globalkey[6], globalkey[7]};
    ak_uint8 m816[8] = {globalkey[8], globalkey[9], globalkey[10], globalkey[11],
                    globalkey[12], globalkey[13], globalkey[14], globalkey[15]};
    ak_uint32 S0 = RSM(m07);
    ak_uint32 S1 = RSM(m816);

    for(int i = 0; i < 2; i++){
        key.Skey[i] = (i==0) ? S0 : S1;
    }
}

void pushKkeys(const ak_uint8 M[16], ak_key_shedule key){

    ak_uint32 Me[2] = {ak_uint8_to_uint32(M[8], M[9], M[10], M[11]),
                  ak_uint8_to_uint32(M[0], M[1], M[2], M[3])};
    ak_uint32 Mo[2] = {ak_uint8_to_uint32(M[12], M[13], M[14], M[15]),
                  ak_uint8_to_uint32(M[4], M[5], M[6], M[7])};
//H
    for(int i=0; i < 20; i++){
        ak_uint32 K2i, K2i1;
        K2i = g(2*i*RHO, Me);
        K2i1 = ROL(g(2*i*RHO + RHO, Mo), 8);
        key.Kkey[2*i] = K2i + K2i1;
        key.Kkey[2*i + 1] = ROL(K2i + 2*K2i1, 9);
    }
}

void ak_create_key_shedule(const ak_uint8 M[16], ak_key_shedule key){
    pushSkeys(&M[16], key);
    pushKkeys(&M[16], key);
}

