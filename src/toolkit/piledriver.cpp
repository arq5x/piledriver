// ***************************************************************************
// bamtools_coverage.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Prints coverage data for a single BAM file 
// ***************************************************************************

#include "piledriver.h"

#include <api/BamMultiReader.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
using namespace BamTools;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
using namespace std;





namespace BamTools {

void profile_coverage(const PileupPosition& pileupData, 
                      vector<SampleCoverage> &sample_cov,
                      map<string, size_t> &sample_map);    

class PileDriverVisitor : public PileupVisitor {
  
    public:
        PileDriverVisitor(const RefVector& references, 
                          ostream* out,
                          int num_samples, 
                          map<string, size_t> &sample_map)
            : PileupVisitor()
            , m_references(references)
            , m_out(out)
            , m_num_samples(num_samples)
            , m_sample_map(sample_map)
        { }
        ~PileDriverVisitor(void) { }
  
    // PileupVisitor interface implementation
    public:
        
        void Header() {
            *m_out << "#chrom\t"
                   << "start\t"
                   << "end\t"
                   << "num_A\t"
                   << "num_C\t"
                   << "num_G\t"
                   << "num_T\t"
                   << "num_D\t"
                   //<< "num_I\t"
                   << "totQ_A\t"
                   << "totQ_C\t"
                   << "totQ_G\t"
                   << "totQ_T\t"
                   << "totQ_D";
                   //<< "totQ_I";
            for (size_t i = 0; i < m_num_samples; ++i)
            {
                *m_out 
                    << "\tsample_" << i + 1;
            }
             *m_out << endl;
        }


        // prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData) {
            
            
            *m_out << m_references[pileupData.RefId].RefName 
                   << "\t" 
                   << pileupData.Position
                   << "\t"
                   << pileupData.Position + 1;
            
            vector<SampleCoverage> sample_cov(m_num_samples + 1);
            
            profile_coverage(pileupData, sample_cov, m_sample_map);
        
            // report the overall coverage for each allele
            // as well as the total quality for each allele
            *m_out << "\t"
                   << sample_cov[m_num_samples].a_cnt << "\t"
                   << sample_cov[m_num_samples].c_cnt << "\t"
                   << sample_cov[m_num_samples].g_cnt << "\t"
                   << sample_cov[m_num_samples].t_cnt << "\t"
                   << sample_cov[m_num_samples].del_cnt << "\t"
                   //<< sample_cov[m_num_samples].ins_cnt << "\t"
                
                   << sample_cov[m_num_samples].a_tot_qual << "\t"
                   << sample_cov[m_num_samples].c_tot_qual << "\t"
                   << sample_cov[m_num_samples].g_tot_qual << "\t"
                   << sample_cov[m_num_samples].t_tot_qual << "\t"
                   << ".";
                   //<< sample_cov[m_num_samples].ins_tot_qual;

            for (size_t i = 0; i < m_num_samples; ++i)
            {
                *m_out 
                  << "\t"
                  // num(A)|totQual(A)
                  << sample_cov[i].a_cnt << "|" << sample_cov[i].a_tot_qual
                  << ","
                  // num(C)|totQual(C)
                  << sample_cov[i].c_cnt << "|" << sample_cov[i].c_tot_qual
                  << ","
                  // num(G)|totQual(G)
                  << sample_cov[i].g_cnt << "|" << sample_cov[i].g_tot_qual
                  << ","
                  // num(T)|totQual(T)
                  << sample_cov[i].t_cnt << "|" << sample_cov[i].t_tot_qual
                  << ","
                  // num(D)|totQual(D)
                  << sample_cov[i].del_cnt << "|.";
                  // << sample_cov[i].del_cnt << "|."
                  // << ","
                  // num(I)|totQual(I)
                  //<< sample_cov[i].ins_cnt << "|" << sample_cov[i].ins_tot_qual;
            }
            *m_out << endl;
        }
        
    private:
        RefVector m_references;
        ostream*  m_out;
        size_t    m_num_samples;
        map<string, size_t> m_sample_map;
};

void profile_coverage(const PileupPosition& pileupData, 
                      vector<SampleCoverage> &sample_cov,
                      map<string, size_t> &sample_map)
{
    size_t all_samples_idx = sample_cov.size() - 1;
    
    for (size_t i = 0; i < pileupData.PileupAlignments.size(); ++i)
    {
        PileupAlignment p_al = pileupData.PileupAlignments[i];
        BamAlignment al = p_al.Alignment;
        
        size_t file_id = sample_map[al.Filename];
        size_t alignmentpos = p_al.PositionInAlignment;
        
        string base;
        if (p_al.IsCurrentDeletion)
            base = "*";
        else    
            base = al.QueryBases[alignmentpos];
        
        short qual =  static_cast<short>(al.Qualities[alignmentpos]) - 33;
        
        if (pileupData.PileupAlignments[i].IsCurrentInsertion)
        {
            sample_cov[file_id].ins_cnt++;
            sample_cov[file_id].ins_tot_qual = qual;
            
            sample_cov[all_samples_idx].ins_cnt++;
            sample_cov[all_samples_idx].ins_tot_qual = qual;
        }
        else {
        
            if ((base == "A") || (base == "a"))
            {
                sample_cov[file_id].a_cnt++;
                sample_cov[file_id].a_tot_qual += qual;
            
                sample_cov[all_samples_idx].a_cnt++;
                sample_cov[all_samples_idx].a_tot_qual += qual;
            }
            else if ((base == "C") || (base == "c"))
            {
                sample_cov[file_id].c_cnt++;
                sample_cov[file_id].c_tot_qual += qual;
            
                sample_cov[all_samples_idx].c_cnt++;
                sample_cov[all_samples_idx].c_tot_qual += qual;
            }
            else if ((base == "G") || (base == "g"))
            {
                sample_cov[file_id].g_cnt++;
                sample_cov[file_id].g_tot_qual += qual;
            
                sample_cov[all_samples_idx].g_cnt++;
                sample_cov[all_samples_idx].g_tot_qual += qual;
            }
            else if ((base == "T") || (base == "t"))
            {
                sample_cov[file_id].t_cnt++;
                sample_cov[file_id].t_tot_qual += qual;
            
                sample_cov[all_samples_idx].t_cnt++;
                sample_cov[all_samples_idx].t_tot_qual += qual;
            }
            else if (base == "*")
            {
                sample_cov[file_id].del_cnt++;
                sample_cov[file_id].del_tot_qual = -1;
            
                sample_cov[all_samples_idx].del_cnt++;
                sample_cov[all_samples_idx].del_tot_qual = -1;
            }
        }
    }
}

} // namespace BamTools

// ---------------------------------------------  
// CoverageSettings implementation

struct PileDriverTool::PileDriverSettings {

    // flags
    bool HasInput;
    bool HasInputFilelist;
    bool HasRegion;

    // filenames
    vector<string> InputFiles;
    string InputFilelist;
    string Region;
    
    // constructor
    PileDriverSettings(void)
        : HasInput(false)
        , HasInputFilelist(false)
        , HasRegion(false)
    { }  

};  

// ---------------------------------------------
// CoverageToolPrivate implementation

struct PileDriverTool::PileDriverToolPrivate {
  
    // ctor & dtor
    public:
        PileDriverToolPrivate(PileDriverTool::PileDriverSettings* settings)
            : m_settings(settings)
            , m_out(cout.rdbuf())
        { }

        ~PileDriverToolPrivate(void) { }
    
    // interface
    public:
        bool Run(void);
        
    // data members
    private: 
        PileDriverTool::PileDriverSettings* m_settings;
        ostream m_out;
        RefVector m_references;
};  

bool PileDriverTool::PileDriverToolPrivate::Run(void) {  
  
    // set to default input if none provided
    if ( !m_settings->HasInput && !m_settings->HasInputFilelist )
        m_settings->InputFiles.push_back(Options::StandardIn());

    // add files in the filelist to the input file list
    if ( m_settings->HasInputFilelist ) {

        ifstream filelist(m_settings->InputFilelist.c_str(), ios::in);
        if ( !filelist.is_open() ) {
            cerr << "bamtools count ERROR: could not open input BAM file list... Aborting." << endl;
            return false;
        }

        string line;
        while ( getline(filelist, line) )
            m_settings->InputFiles.push_back(line);
    }

    // open reader without index
    BamMultiReader reader;
    if ( !reader.Open(m_settings->InputFiles) ) {
        cerr << "bamtools count ERROR: could not open input BAM file(s)... Aborting." << endl;
        return false;
    }

    // retrieve references
    m_references = reader.GetReferenceData();
    
    // set up our output 'visitor'
    map<string, size_t> sample_map;
    for (size_t i = 0; i < m_settings->InputFiles.size(); ++i)
        sample_map[m_settings->InputFiles[i]] = i;
    
    PileDriverVisitor* cv = 
        new PileDriverVisitor(m_references, 
                              &m_out,
                              m_settings->InputFiles.size(),
                              sample_map);
    
    // print a header
    cv->Header();
    
    // set up pileup engine with 'visitor'
    PileupEngine pileup;
    pileup.AddVisitor(cv);
    
    // process input data
    BamAlignment al;    
    while ( reader.GetNextAlignment(al) )
    {  
        pileup.AddAlignment(al);
    }
    // clean up 
    reader.Close();
    //if ( m_settings->HasOutputFile )
    //    outFile.close();
    delete cv;
    cv = 0;
    
    // return success
    return true;
}

// ---------------------------------------------
// CoverageTool implementation
/*
PileDriverTool::PileDriverTool(void) 
    : AbstractTool()
    , m_settings(new CoverageSettings)
    , m_impl(0)
{ 
    // set program details
    Options::SetProgramInfo("bamtools piledriver", "prints pileup data for a single BAM file", "[-in <filename>] [-out <filename>]");
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",  "BAM filename", "the input BAM file", "", m_settings->HasInputFile,  m_settings->InputBamFilename, IO_Opts, Options::StandardIn());
    Options::AddValueOption("-out", "filename",     "the output file",    "", m_settings->HasOutputFile, m_settings->OutputFilename,   IO_Opts, Options::StandardOut());
}
*/



PileDriverTool::PileDriverTool(void) 
    : AbstractTool()
    , m_settings(new PileDriverSettings)
    , m_impl(0)
{ 
    // set program details
    Options::SetProgramInfo("bamtools piledriver", "prints number of alignments in BAM file(s)",
                            "[-in <filename> -in <filename> ... | -list <filelist>] [-region <REGION>]");
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",     "BAM filename", "the input BAM file(s)", "", m_settings->HasInput,  m_settings->InputFiles, IO_Opts, Options::StandardIn());
    Options::AddValueOption("-list",   "filename", "the input BAM file list, one line per file", "", m_settings->HasInputFilelist,  m_settings->InputFilelist, IO_Opts);
    Options::AddValueOption("-region", "REGION",
                            "genomic region. Index file is recommended for better performance, and is used automatically if it exists. See \'bamtools help index\' for more details on creating one",
                            "", m_settings->HasRegion, m_settings->Region, IO_Opts);
}



PileDriverTool::~PileDriverTool(void) { 

    delete m_settings;
    m_settings = 0;
    
    delete m_impl;
    m_impl = 0;
}

int PileDriverTool::Help(void) { 
    Options::DisplayHelp();
    return 0;
} 

int PileDriverTool::Run(int argc, char* argv[]) { 

    // parse command line arguments
    Options::Parse(argc, argv, 1);
    
    // initialize CoverageTool with settings
    m_impl = new PileDriverToolPrivate(m_settings);
    
    // run CoverageTool, return success/fail
    if ( m_impl->Run() ) 
        return 0;
    else 
        return 1;
}  
