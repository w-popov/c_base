#include "gF8.hpp"
#include "HW9.h"
#include <sstream>
#include <algorithm>
#include <FL/Fl_Button.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Grid.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl.H>

#ifdef _WIN32
void seticon (Fl_Window *win)
{
    HWND hwnd = fl_xid(win);
    if (hwnd)
    {
        HICON hIcon = LoadIcon(GetModuleHandle(NULL), "MAINICON");
        if (hIcon)
        {
            SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hIcon);
            SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hIcon);
        }
    }
}
#endif

inline F8::F8(const char *arr_values)
    : array(std::move(parse_input(arr_values))){};

F8::F8(F8 &&other) noexcept : array(std::move(other.array)){};

F8::F8(const F8 &other) : array(other.array){};

std::vector<int> &F8::init(const char *input)
{
    auto parsed = parse_input(input);
    array = std::move(parsed);
    return array;
}

bool F8::is_empty() const
{
    return array.empty();
}

std::vector<int> F8::parse_input(const char *input)
{
    std::vector<int> v;
    if (input == nullptr || *input == '\0')
    {
        return v;
    }
    std::istringstream iss(input);
    int number;
    while (iss >> number && number)
    {
        v.push_back(number);
    }
    return v;
}

int F8::missing_num()
{
    return missing_number(array.data(), static_cast<int>(array.size()));
}

const std::vector<int> &F8::get_array() const
{
    return array;
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

WindowF8::WindowF8(int w, int h, const char *title)
    : Fl_Double_Window(w, h, title)
{
    Fl_Grid *main_grid = new Fl_Grid(0, 0, w, h);
    main_grid->layout(2, 1, 5, 5); // 2 строки, 1 колонка
    main_grid->color(fl_rgb_color(250, 250, 250));

    /* ВЕРХНЯЯ ЧАСТЬ (Кнопки и поля ввода) */
    Fl_Grid *top_grid = new Fl_Grid(0, 0, w, 290);
    top_grid->layout(6, 6, 5, 5);
    top_grid->color(fl_rgb_color(250, 250, 250));

    info_box = new Fl_Box(0, 0, w, 80, WindowF8::info);
    info_box->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE | FL_ALIGN_WRAP);
    info_box->labelfont(FL_HELVETICA);
    info_box->labelsize(13);
    info_box->color(FL_FREE_COLOR);
    info_box->box(FL_NO_BOX);
    top_grid->widget(info_box, 0, 0, 1, 6, FL_GRID_FILL);
    top_grid->row_gap(0, 5);

    in = new Fl_Input(0, 0, 550, 30);
    in->textsize(14);
    top_grid->widget(in, 2, 2, 1, 3, FL_GRID_HORIZONTAL);
    top_grid->row_gap(2, 10);

    out = new Fl_Output(0, 0, 550, 30);
    out->color(FL_LIGHT1);
    out->textsize(14);
    top_grid->widget(out, 3, 2, 1, 3, FL_GRID_HORIZONTAL);
    top_grid->row_gap(3, 10);

    btn = new Fl_Button(0, 0, 200, 40, "Сортировать");
    btn->callback(button_callback, this);
    top_grid->widget(btn, 4, 2, 1, 1, FL_GRID_HORIZONTAL);
    top_grid->col_gap(2, 20);

    btn_clear = new Fl_Button(0, 0, 200, 40, "Очистить");
    btn_clear->callback(button_clear, this);
    top_grid->widget(btn_clear, 4, 4, 1, 1, FL_GRID_HORIZONTAL);
    top_grid->col_gap(2, 20);

    int cw[] = {20, 0, 0, 20, 20, 20};
    int rw[] = {20, 0, 0, 0, 0, 10};
    top_grid->col_weight(cw, 6);
    top_grid->row_weight(rw, 6);
    top_grid->end();

    /* НИЖНЯЯ ЧАСТЬ (Сетка для ячеек) */
    scroll = new Fl_Scroll(0, 0, w, h - 290);
    scroll->scrollbar.size(30, scroll->scrollbar.h());
    scroll->scrollbar.align(FL_ALIGN_RIGHT);
    scroll->color(fl_rgb_color(250, 250, 250));
    scroll->when(FL_WHEN_NEVER);
    cells_grid = new Fl_Grid(0, 0, w - 20, 100); 
    cells_grid->layout(1, 1, 5, 5);
    cells_grid->color(fl_rgb_color(250, 250, 250));
    cells_grid->end();

    scroll->begin();
    scroll->add(cells_grid);
    scroll->end();
    scroll->resizable(cells_grid);

    main_grid->begin();
    main_grid->widget(top_grid, 0, 0, 1, 1, FL_GRID_FILL);
    main_grid->widget(scroll, 1, 0, 1, 1, FL_GRID_FILL);
    main_grid->end();

    this->end();
    this->resizable(main_grid);
    this->size_range(800, 600, 900, 800);
}

void WindowF8::print_output()
{
    int missing = f8.missing_num();
    out->value(std::to_string(missing).c_str());
    create_cells();
}

void WindowF8::button_callback(Fl_Widget *, void *data)
{
    WindowF8 *win = static_cast<WindowF8 *>(data);
    win->f8.init(win->in->value());
    if (win->f8.is_empty())
    {
        win->out->value("Введите числа для сортировки.");
        return;
    }
    win->print_output();
}

void WindowF8::button_clear(Fl_Widget *, void *data)
{
    WindowF8 *win = static_cast<WindowF8 *>(data);
    win->in->value("");
    win->out->value("");
    win->clear_cells();
}

void WindowF8::clear_cells()
{
    if (cells_grid)
    {
        cells_grid->begin();
        for (auto *cell : cell_widgets)
        {
            cells_grid->remove(cell);
            delete cell;
        }
        cell_widgets.clear();
        cells_grid->end();
        cells_grid->size(cells_grid->w(), 100);
        cells_grid->redraw();
        scroll->redraw();
    }
}

void WindowF8::create_cells()
{
    clear_cells();

    const auto &arr = f8.get_array();
    if (arr.empty())
    {
        return;
    }

    int missing = f8.missing_num();
    std::vector<int> sorted = arr;
    std::sort(sorted.begin(), sorted.end());

    int cols = 5;
    int cell_w = 40;
    int cell_h = 40;
    int gap = 5;

    int total_items = sorted.size() + (missing > 0 ? 1 : 0);
    int rows = (total_items + cols - 1) / cols;

    cells_grid->layout(rows, cols, gap, gap);
    cells_grid->begin();

    int idx = 0;

    /* Создать ячейку пропущенного числа (если она есть) */
    if (missing > 0)
    {
        int r = idx / cols; // r = 0
        int c = idx % cols; // c = 0

        Fl_Box *missing_cell = new Fl_Box(0, 0, cell_w, cell_h);
        missing_cell->size(cell_w, cell_h);

        char buf[16];
        snprintf(buf, sizeof(buf), "%d", missing);
        missing_cell->copy_label(buf);
        missing_cell->box(FL_BORDER_BOX);
        missing_cell->color(fl_rgb_color(30, 161, 26));
        missing_cell->labelfont(FL_COURIER_BOLD);
        missing_cell->labelsize(15);
        missing_cell->labelcolor(FL_WHITE);

        cells_grid->widget(missing_cell, r, c, 1, 1, FL_GRID_FILL);
        cell_widgets.push_back(missing_cell);
        /* Увеличить счётчик, чтобы следующие ячейки пошли со 2-й позиции */
        idx++; 
    }

    /* Создать все остальные отсортированные ячейки */
    for (int val : sorted)
    {
        int r = idx / cols;
        int c = idx % cols;

        Fl_Box *cell = new Fl_Box(0, 0, cell_w, cell_h);
        cell->size(cell_w, cell_h);

        char buf[16];
        snprintf(buf, sizeof(buf), "%d", val);
        cell->copy_label(buf);
        cell->box(FL_BORDER_BOX);
        cell->labelfont(FL_COURIER_BOLD);
        cell->labelsize(14);

        if (val == missing - 1)
        {
            cell->color(fl_rgb_color(255, 200, 200));
        }
        else if (val == missing + 1)
        {
            cell->color(fl_rgb_color(200, 255, 200));
        }
        else
        {
            cell->color(FL_WHITE);
        }

        cells_grid->widget(cell, r, c, 1, 1, FL_GRID_FILL);
        cell_widgets.push_back(cell);
        idx++;
    }

    cells_grid->end();

    int new_height = rows * (cell_h + gap) + gap;
    cells_grid->size(cells_grid->w(), new_height);

    scroll->redraw();
    cells_grid->redraw();
}

int main (int argc, char **argv)
{
    Fl::scheme("gtk+");
    WindowF8 *win = new WindowF8(800, 700, "ДЗ-9: F8");
    win->show(argc, argv);

    Fl::check();

#ifdef _WIN32
    seticon(win);
#endif
    return Fl::run();
}
