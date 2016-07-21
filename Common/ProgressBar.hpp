/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/** \file ProgressBar.hpp
 *  \brief Write Me!!!!!!
 *  \author Ryan M. Richard
 *  \version 1.0
 *  \date July 18, 2016
 */

#ifndef PULSAR_GUARD_PROGRESSBAR_HPP
#define PULSAR_GUARD_PROGRESSBAR_HPP

#include<iostream>

/** \brief Prints a progress bar to show progress during long, ardourous
 *         tasks
 *
 *   The bar for a task that is between 42% and 44% complete would look like:
 *   \verbatim
 *   0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%
 *   |----|----|----|----|----|----|----|----|----|----|
 *   **********************
 *   \endverbatim
 *
 *   We need to print 50 characters. If we have NTasks%50==0 we simply
 *   print a character every time we've incremented Increment_ times. If
 *   NTasks%5=R, the first R times we increment to Increment_ we actually
 *   make it increment to Increment_+1.  This means the first R increments
 *   will be slower than the latter NTasks-R increments.  My philosophy is
 *   you become more impatient as you near the end, so make those run faster.
 *
 *   If we have less than 50 tasks we have to print multiple characters
 *   at a time.  In this case we print floor(50/NTasks) characters every
 *   call for the first NTask-50%NTasks calls, and floor(50/NTasks)+1 characters
 *   on the last call.  Again, we speed-up at the end.
 *
 *   The progress bar will get garbled in the output if the thing you
 *   are doing also prints to the output file.  The progress bar does
 *   not place a newline character even after finishing, thus you need
 *   to do it.  This prevents the output from being garbled by your
 *   progress bar not being called enough times (possible if you're
 *   estimating the number of iterations)
 *
 */
class ProgressBar{
   private:
      typedef size_t NTask_t;
      ///Flag telling us whether we are above or below 50 Tasks
      bool LessThan50_;
      ///The number of increments that have printed so far
      NTask_t NStars_;
      ///The number of tasks
      NTask_t NTasks_;
      ///The number of tasks before we increment
      NTask_t Increment_;
      ///The non-evenly distributed
      NTask_t Remainder_;
      ///Our current Task
      NTask_t Current_;
      ///The number of characters to print at a time
      int NChars_;
      ///The character we are printing
      char Char_;
      ///The ostream we are printing to
      std::ostream& os_;
   public:
      ///Makes a progress bar and prints the non-bar part
      ProgressBar(const NTask_t NTasks,std::ostream& os);
      ///Increment and print if we've done enough increments
      ProgressBar& operator++();
};

#endif /* PULSAR_GHUARD_PROGRESSBAR_HPP */

