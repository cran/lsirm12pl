// [[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <cstdio>
#include <chrono>
#include <Rcpp.h>

using namespace Rcpp;

void progressbar(int step, int total)
{
  // progress width
  const int pwidth = 72;

  // minus label len
  int pos = (step * pwidth) / total;
  int percent = (step * 100) / total;

  // calculate elapsed time in seconds
  static auto start_time = std::chrono::steady_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count();

  // calculate remaining time in seconds
  auto remaining_seconds = elapsed_seconds * (total - step) / step;
  int hours = remaining_seconds / 3600;
  int minutes = (remaining_seconds % 3600) / 60;
  int seconds = remaining_seconds % 60;

  // calculate total elapsed time in seconds
  // int total_elapsed_seconds = elapsed_seconds + remaining_seconds;
  // int total_hours = total_elapsed_seconds / 3600;
  // int total_minutes = (total_elapsed_seconds % 3600) / 60;
  // int total_seconds = total_elapsed_seconds % 60;

  // fill progress bar with =
  Rcpp::Rcout << "[";
  for(int i = 0; i < pos; i++){
    //std::printf("%c", '=');
    Rprintf("=");
  }

  // fill progress bar with spaces
  Rprintf("%*c", pwidth - pos + 1, ']');

  // print percentage, ETA, and total elapsed time
  Rprintf(" %3d%% ETA: %02d:%02d:%02d \r", percent, hours, minutes, seconds);

  // flush output to make sure it's displayed immediately
  // std::cout.flush();
  R_FlushConsole();
}
