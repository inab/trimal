//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_READWRITEMANAGER_H
#define TRIMAL_READWRITEMANAGER_H
#include <ostream>
using namespace std;

class newAlignment;

class ReadWriteManager {

public:
    int getInputFormat(void);

    int getOutputFormat(void);

    bool loadAlignment(char *alignmentFile);

    bool saveAlignment(char *destFile);

    bool printAlignment(void);

    int formatInputAlignment(char *);

    int formatInputFile(void);

    int typeInputFile(void);

    bool loadPhylipAlignment(char *);

    bool loadFastaAlignment(char *);

    bool loadClustalAlignment(char *);

    bool loadNexusAlignment(char *);

    bool loadMegaInterleavedAlignment(char *);

    bool loadMegaNonInterleavedAlignment(char *);

    bool loadNBRF_PirAlignment(char *);

    bool loadPhylip3_2Alignment(char *);

    void alignmentClustalToFile(ostream &);

    void alignmentNBRF_PirToFile(ostream &);

    void alignmentFastaToFile(ostream &);

    void alignmentPhylip3_2ToFile(ostream &);

    void alignmentPhylipToFile(ostream &);

    void alignmentPhylip_PamlToFile(ostream &);

    void alignmentNexusToFile(ostream &);

    void alignmentMegaToFile(ostream &);

    bool alignmentSummaryHTML(char *, int, int, int *, int *, float *);

    bool alignmentColourHTML(ostream &);

    void getSequences(ostream &);

private:

    friend class newAlignment;

    newAlignment* _alignment;

    ReadWriteManager(newAlignment* parent);
};


#endif //TRIMAL_READWRITEMANAGER_H
