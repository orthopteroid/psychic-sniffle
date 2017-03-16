// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license

#ifndef PROJECT_CPUINFO_H
#define PROJECT_CPUINFO_H

// http://stackoverflow.com/a/3082553/968363
// but dividing cores by 2 when HTT is active
static inline int EnumCores()
{
    const int dwIntel = 'uneG'; // GenuineIntel
    const int dwAMD = 'htuA'; // AuthenticAMD

    auto cpuID = [](unsigned i, unsigned regs[4]) -> void {
        asm volatile
        ("cpuid" : "=a" (regs[0]), "=b" (regs[1]), "=c" (regs[2]), "=d" (regs[3])
        : "a" (i), "c" (0));
    };

    unsigned regs[4];

    cpuID(1, regs);
    unsigned cpuFeatures = regs[3]; // EDX

    cpuID(1, regs);
    unsigned logical = (regs[1] >> 16) & 0xff; // EBX[23:16]
    unsigned cores = logical;

    cpuID(0, regs);
    if (regs[1] == dwIntel) {
        cpuID(4, regs);
        cores = ((regs[0] >> 26) & 0x3f) + 1; // EAX[31:26] + 1

    } else if (regs[1] == dwAMD) {
        cpuID(0x80000008, regs);
        cores = ((unsigned)(regs[2] & 0xff)) + 1; // ECX[7:0] + 1
    }

    bool hyperThreads = (cpuFeatures & (1 << 28)) && cores < logical;

    // works for me: an older i7 920
    return cores / (hyperThreads ? 2 : 1);
}

#endif //PROJECT_CPUINFO_H
