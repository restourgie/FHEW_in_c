#include "params.h"

const ZmodQ vgprime[3] = {v, 2048, 4194304};
const int g_bits[3] = {11, 11, 10};
const int g_bits_32[3] = {21, 21, 22};

const ZmodQ KS_table[7] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25,
			   25*25*25*25*25*25};
const int BS_table[2] = {1,23};