// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license
//
// Splices two byte arrays of the same length together at the specified bit.
// When splicing, bits from a will be written to lower memory than bits from b.
// In the transition byte, high bits come from a and high bits come from b. Does it matter?

#ifndef PROJECT_SPLICE_H
#define PROJECT_SPLICE_H

namespace util {

template<uint Size>
void splice(uint8_t *out, uint8_t *a, uint8_t *b, uint uRand) {
    const uint8_t u8Mask[] = {
        0b00000000, 0b00000001, 0b00000011, 0b00000111,
        0b00001111, 0b00011111, 0b00111111, 0b01111111
    };

    uint uBit = uRand % ( Size * 8 );
    int i = 0;
    while (i < (uBit / 8)) out[i++] = a[i];
    out[i++] = (a[i] & ~u8Mask[uBit & 7]) | (b[i] & u8Mask[uBit & 7]);
    while (i < Size) out[i++] = b[i];
}

}

#endif //PROJECT_SPLICE_H
