#include "pulsar_modules/common/ProgressBar.hpp"
#include <pulsar/math/Cast.hpp>
#include <cmath>

ProgressBar::ProgressBar(size_type NTasks,std::ostream& os):
      LessThan50_(NTasks<50),NStars_(0),NTasks_(NTasks),Current_(0),
       NChars_(1),Char_('*'),os_(os){
   if(!LessThan50_){
      Remainder_=(NTasks%50);
      Increment_=(NTasks-Remainder_)/50;
   }
   else{
      Remainder_=50%NTasks;
      NChars_=pulsar::numeric_cast<size_type>(std::floor(50/NTasks));
      Increment_=1;
   }
   os_<<"0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%"<<std::endl;
   os_<<"|----|----|----|----|----|----|----|----|----|----|"<<std::endl;
   os_<<Char_;
   os_.flush();
}

ProgressBar& ProgressBar::operator++(){
   ++Current_;
   if(LessThan50_){
      os_<<std::string(NChars_+
            (NStars_>NTasks_-Remainder_-1?1:0),Char_);
      NStars_++;
   }
   else if((Current_==Increment_+1&&Remainder_>0)||
      (Current_==Increment_&&Remainder_==0)){
      os_<<std::string(NChars_,Char_);
      if(Remainder_>0)Remainder_--;
      Current_=0;
      NStars_++;
   }
   os_.flush();
   //else the bar has been printed
   return *this;
}
