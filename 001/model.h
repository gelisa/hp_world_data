#ifndef __SPECIE_H
#define __SPECIE_H

#include <vector>
#include <string>
#include <list>
#include "parameter.h"
#include "reaction.h"

class Specie
{
public:
    std::string modelName;
    Specie(){};
    Specie(std::string id);
    ~Specie();
    
    //attributes
    std::string m_id;
    int m_length;
    float m_hydrophobicity;
    
    //constants
    
    
    //methods
    std::list<Reaction> reactions(Specie specie);
    
    friend std::ostream& operator<<(std::ostream& os, const Specie& sp);
    
private:
    //imports either H* or P* 
    void importHorP( std::list<Reaction>& allReactions,Specie other,
                       float impRate,std::string HorP);

    void degradeIt(std::list<Reaction>& allReactions,Specie specie,
              float degrRate);

    void aggregateIt(std::list<Reaction>& allReactions,Specie specie,
                     float aggRate,float aggPower);
    void hydrolyseIt(std::list<Reaction>& allReactions,Specie specie,
                     float dH);
    void growIt(std::list<Reaction>& allReactions,Specie specie,
                float alpha, int maxLength);
    void growOther(std::list<Reaction>& allReactions,Specie specie,
                float alpha, int maxLength);
    
};


#endif // __SPECIE_H
