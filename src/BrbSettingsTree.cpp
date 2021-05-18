
#include "BrbSettingsTree.hxx"

BrbSettingsTree* BrbSettingsTree::singleton_= nullptr;;

/**                                                                                                                                                                    
 * Static methods should be defined outside the class.                                                                                                                 
 */
BrbSettingsTree *BrbSettingsTree::Get()
{

  if(singleton_==nullptr){
    singleton_ = new BrbSettingsTree();
  }
  return singleton_;
}



int BrbSettingsTree::LoadSettingsTree(TFile *file){

  if(!file){
    std::cerr << "LoadSettingsTree; can't open file "  << std::endl;
    return -1;
  }
  
  TTree* settings_tree = (TTree*)file->Get("settings_tree;1");
  
  if(!settings_tree){
    std::cerr << "Failed to load settings_tree from file "<< std::endl;
    file->Close();
    return -1;
  }
  std::cout << "Opened BRB settings_tree from file " << std::endl;

  // Get baselines from tree
  TBranch* brBaseline = settings_tree->GetBranch("CalcBaseline");
  if(!brBaseline){
    std::cerr << "Can't read branch CalcBaseline" << std::endl;
    return -1;
  }
  double thisbaseline[20];
  brBaseline->SetAddress(&thisbaseline);

  // Get HVSetpoints from tree
  TBranch* brHV = settings_tree->GetBranch("HVsetpoints");
  if(!brHV){
    std::cerr << "Can't read branch HVSetpoints" << std::endl;
    return -1;
  }
  double ihv[20];
  brHV->SetAddress(&ihv);

  // Get the single event in tree
  settings_tree->GetEvent(0);

  // Store values 
  for(int i = 0; i < 20 ; i++){

    // convert baseline from counts to voltage
    fBaselines.push_back(thisbaseline[i] * (2.0/4096.0));

    fHV_setpoints.push_back(ihv[i]);
    
  }

  
  return 0;

}
