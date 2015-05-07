#ifndef COMMON_H
#define COMMON_H

#include "../LWE.h"
#include "../FHEW.h"


void SaveSecretKey(SecretKey* ct, char* filepath);
SecretKey* LoadSecretKey(char* filepath);


void SaveEvalKey(const EvalKey *EK, char* filepath);
EvalKey* LoadEvalKey(char* filepath);


void SaveCipherText(const CipherText* ct, char* filepath);
CipherText* LoadCipherText(char* filepath);


#endif