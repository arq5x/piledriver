// ***************************************************************************
// bamtools_convert.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 July 2010
// ---------------------------------------------------------------------------
// Converts between BAM and a number of other formats
// ***************************************************************************

#ifndef PILE_H
#define PILE_H

#include <vector>
#include "bamtools_tool.h"


namespace BamTools { 
  
class PileDriverTool : public AbstractTool {
  
    public:
        PileDriverTool(void);
        ~PileDriverTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private: 
        struct PileDriverSettings;
        PileDriverSettings* m_settings;
        
        struct PileDriverToolPrivate;
        PileDriverToolPrivate* m_impl;
};

struct SampleCoverage {
    size_t a_cnt;
    size_t c_cnt;
    size_t g_cnt;
    size_t t_cnt;
    size_t del_cnt;
    size_t ins_cnt;

    size_t a_tot_qual;
    size_t c_tot_qual;
    size_t g_tot_qual;
    size_t t_tot_qual;
    int del_tot_qual;
    std::vector<std::string> ins_alleles;
};
  
} // namespace BamTools

#endif // PILE_H
