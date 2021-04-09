cmake_minimum_required(VERSION 2.8) # установка минимальной версии cmake
project(FFT) # название проекта

set(SOURCES main.cpp DFT.cpp) # Создание переменной
add_executable(FFT ${SOURCES}) # Объявляет исполняемый модуль с именем sample_app
