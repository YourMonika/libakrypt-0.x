/* ----------------------------------------------------------------------------------------------- */
/*                 Тестирование алгоритма блочного шифрования TwoFish (128-bit)                    */
/* ----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <libakrypt.h>

int main(){
    ak_uint32 plain[4] = {0, 0, 0, 0};
    ak_uint8 M[16] = {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0};
    /*
     * Mute BYTE-MASSIVE M to UINT-MASSIVE
     */
    ak_uint32 gM[4] = {ak_uint8_to_uint32(M[0], M[1], M[2], M[3]),
                       ak_uint8_to_uint32(M[4], M[5], M[6], M[7]),
                       ak_uint8_to_uint32(M[8], M[9], M[10], M[11]),
                       ak_uint8_to_uint32(M[12], M[13], M[14], M[15])};
    ak_key_shedule key_shedule;

    printf("Creating Key Shedule \n");
    ak_create_key_shedule((const ak_uint8 *) M, key_shedule);
    for(int j = 0; j < 40; j++){
        printf("Key Shedule [%i] K-key = %lu\n", j, key_shedule.Kkey[j]);
    }
    printf("Finished Key SHedule creating, start encryption\n");

    ak_uint32 chipher[4];
    ak_twofish_encrypt(plain, key_shedule, chipher);
    printf("CHIPHER: \n");
    for(int j = 0; j < 4; j++){
        printf("%i pos - %i \n", j, (ak_uint32) chipher[j]);
    }

    return 0
}