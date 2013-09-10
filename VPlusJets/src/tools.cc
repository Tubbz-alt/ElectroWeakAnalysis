/*
 * =====================================================================================
 * 
 *       Filename:  tools.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08/07/13 02:04:15 CDT
 *       Revision:  none
 *       Compiler:  gcc, root
 * 
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 * 
 * =====================================================================================
 */

#include "ElectroWeakAnalysis/VPlusJets/interface/tools.h"

void print_p4(fastjet::PseudoJet tmpJ, std::string tmpName,bool extra_info){
	if(extra_info){ 
		std::cout<<tmpName<<" PseudoJet(pt,eta,phi,E,m,pdgID,charge)=("<<tmpJ.pt()<<","<<tmpJ.eta()<<","<<tmpJ.phi()<<","<<tmpJ.E()<<","<<tmpJ.m()<<","<<tmpJ.user_info<PseudoJetUserInfo>().pdg_id()<<","<<tmpJ.user_info<PseudoJetUserInfo>().charge()<<")"<<std::endl;
	}else{
		std::cout<<tmpName<<" PseudoJet(pt,eta,phi,E,m)=("<<tmpJ.pt()<<","<<tmpJ.eta()<<","<<tmpJ.phi()<<","<<tmpJ.E()<<","<<tmpJ.m()<<")"<<std::endl;
	}
}

void BREAK(std::string info){ 
	std::cout<<info<<std::endl<<"Enter a char to continue..."<<std::endl;
	char tmp;std::cin>>tmp;
}

