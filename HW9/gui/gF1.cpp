#include "gF1.hpp"
#include <string>
#include <FL/Fl.H>
#include <FL/fl_draw.H>

MyWindow::MyWindow(int w, int h, const char *title)
    : Fl_Double_Window(w, h, title), counter(0)
{
    begin();
    btn = new Fl_Button(100, 50, 150, 40, "Button нажми");
    btn->labelsize(14); 
    btn->callback(btn_cb, this);
    out = new Fl_Output(100, 120, 150, 30);
    out->value("0");
    out->textsize(14);
    
    end();
}

void MyWindow::btn_cb(Fl_Widget *, void *data)
{
    MyWindow *win = static_cast<MyWindow *>(data);
    win->counter++;
    win->out->value(std::to_string(win->counter).c_str());
}

int main ()
{
    Fl::set_font(FL_HELVETICA, "DejaVu Sans");
    MyWindow win(350, 200, "ДЗ-9. с GUI");
    win.show();
    return Fl::run();
}

