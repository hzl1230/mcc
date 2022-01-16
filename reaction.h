#ifndef _REACTION_ 
#define _REACTION_ 

#include "parse.h"
#include "espic_info.h"
#include "espic_type.h"
#include "espic_math.h"
#include <algorithm>
#include <iostream>

typedef std::vector<std::string> StringList;
typedef std::pair<std::string,std::string> ReactPair;

class Reaction
{
public:
    Reaction (std::string file, const ReactPair& spair, int id);
    ~Reaction ();

    std::vector<Real> en_cs(Real en)
    {
        // Get cs_info through 1DLinearInterpoltationMethod
        if(en <= de_) { 
            info ={0.,0.,0.};
            return info;
        } 

        Real deinv(1./de_);
        Real ei = en*deinv - 1;
        
        // int elo = static_cast<int>(ei);
        int elo = static_cast<int>(ei);
        
        if (elo > arr_length) {
            std::cout << en << " " << elo << std::endl;
            std::string var_name("info_arr_size");
            espic_error(out_bound_info(var_name, infile));
        }
        Real wgt[2], temp;
        // Cal node weight 
        wgt[1] = ei - elo;
        wgt[0] = 1 - wgt[1];
        
        if(!info.empty()) info.clear();
        for(int i = 0; i < info_size; ++i) {
            temp = wgt[0]*info_arr[elo][i] + wgt[1]*info_arr[elo+1][i];
            info.emplace_back(std::move(temp));
        }
        return info;
    }
    void find_max_coll_freq();

    std::vector<std::vector<int> > prodid_arr;
    bool is_background_collision;

    const int size() const { return arr_length; }
    const int isize() const { return info_size; }
    const std::vector<Real>& th() const { return threshold; }
    const ReactPair& pair() const { return spec_pair; }
    const int r_index() { return reaction_id; }
    const Real de() { return de_; }
    const std::vector<Real>& csection(int i) { return info_arr[i]; }
    const StringList& get_types() { return types; }
    
    const std::string get_file() const { return infile; }
    const StringList& get_prod(int i) const { return product_arr[i]; }

    const int sub_cycle() { return n_sub; } 
    Real& mr() { return mr_; }
    const Real& mr() const { return mr_; }
    Real& max_coll_freq() { return nu_max; }
    const Real& max_coll_freq() const { return nu_max; }
    
private:
    std::string infile;
    int info_size;
    const ReactPair spec_pair;
    int reaction_id, arr_length;
    Real de_, deinv_;
    std::vector<Real> threshold;
    std::vector<Real> info;
    std::vector<std::vector<Real>> info_arr; // info_arr of every energy
    std::vector<Real> energy;                // energy_arr of total energy bin
    StringList types;
    std::vector<StringList> product_arr;
    int n_sub;
    Real mr_;
    Real nu_max;

    void element_resize();
    void resize_threshold();

};

inline Real toReal(const std::string& str) 
{
    return (Real)atof(str.c_str());
}

#endif
