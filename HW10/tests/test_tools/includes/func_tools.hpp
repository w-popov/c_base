#pragma once

template <typename T>
int* fill_array (T& array)
{
    for (auto &item : array) {
        if (scanf("%d", &item) != 1)
        {
            printf("Error scanf in fill_array()\n");
            exit(EXIT_FAILURE);
        }
    }
    return array.data();
}
