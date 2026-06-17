#pragma once

#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Output.H>

class MyWindow : public Fl_Double_Window {
    Fl_Button *btn;
    Fl_Output *out;
    int counter;

    static void btn_cb(Fl_Widget *w, void *data);
public:
    MyWindow(int w, int h, const char *title = "ДЗ-9. FLTK GUI");
};