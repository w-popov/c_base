#pragma once

#include <FL/Fl_Double_Window.H>
#include <vector>

class Fl_Multiline_Output;
class Fl_Input;
class Fl_Output;
class Fl_Button;

constexpr const char *INFO = "ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г\n\
F1: Написать только одну функцию, которая сортирует массив по возрастанию void sort_array(int size, int a[])\n";

class MyWindow : public Fl_Double_Window
{
  private:
    Fl_Multiline_Output *info;
    Fl_Input *in;
    Fl_Output *out;
    Fl_Button *btn;
    std::vector<int> array;

    static void btn_cb (Fl_Widget *w, void *data);
    std::vector<int>& parse_input(const char *input);

  public:
    MyWindow(int w, int h, const char *title = "ДЗ-9: F1");

    void sort_arr (int *array, int size);
};