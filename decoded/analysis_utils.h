#pragma once

/*******************************************************************************/

struct data_t {
  int fifo;
  int type;
  int counter;
  int column;
  int pixel;
  int tdc;
  int rollover;
  int coarse;
  int fine;
};

enum type_t {
  kAlcorHit = 1,
  kTriggerTag = 9,
  kStartSpill = 7,
  kEndSpill = 15
};

/*******************************************************************************/

std::pair<int, int> get_index(int pixel, int column);
int get_eochannel(int pixel, int column);

/*******************************************************************************/

int eo2do[32] = {22, 20, 18, 16, 24, 26, 28, 30, 25, 27, 29, 31, 23, 21, 19, 17,
                 9, 11, 13, 15, 7, 5, 3, 1, 6, 4, 2, 0, 8, 10, 12, 14};

int
get_dochannel(int pixel, int column)
{
  int eoch = pixel + 4 * column;
  int doch = eo2do[eoch];
  return doch;
}

/*******************************************************************************/

std::pair<int, int>
get_index(int pixel, int column)
{
  int doch = get_dochannel(pixel, column);
  int ix = doch / 4;
  int iy = doch % 4;
  return {ix, iy};
}

