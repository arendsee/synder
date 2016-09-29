#ifndef __ARGUMENTS_H__
#define __ARGUMENTS_H__

typedef enum {
    C_UNSET,
    C_FILTER,
    C_DEBUG,
    C_DUMP,
    C_MAP,
    C_COUNT,
    C_SEARCH
} Command;

class Arguments
{
public:
    FILE *synfile  = nullptr;
    FILE *intfile  = nullptr;
    FILE *tclfile  = nullptr;
    FILE *qclfile  = nullptr;
    bool swap      = false;
    long k         = 0;
    char trans     = 'i';

    ~Arguments()
    {
        if (synfile != nullptr)
            fclose(synfile);
        if (intfile != nullptr)
            fclose(intfile);
        if (tclfile != nullptr)
            fclose(tclfile);
        if (qclfile != nullptr)
            fclose(qclfile);
    }
};

#endif
