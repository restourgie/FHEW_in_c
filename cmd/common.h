#ifndef COMMON_H
#define COMMON_H

#include "../LWE.h"
#include "../FHEW.h"


void SaveSecretKey(SecretKey LWEsk, char* filepath);
SecretKey* LoadSecretKey(char* filepath);


void SaveEvalKey(EvalKey *EK, char* filepath);
EvalKey* LoadEvalKey(char* filepath);


void SaveCipherText(CipherText* ct, char* filepath);
CipherText* LoadCipherText(char* filepath);


#endif