# DNA Motif Finder

Вычисление количества мотивов в наборе цепочек ДНК с использованием MPI+OpenMP.

## Постановка задачи

1. В заданном входном текстовом файле содержатся массивы фиксированной длины (40 элементов) - ChIP-seq, состоящие из кодов аминокислот ДНК - буквы A,T,G,C. Пример файла находится в материалах - "FoxA2_5000.fst".

2. Задан  мотив - шаблон из 8 позиций. Отличается от обычного массива тем, что на позиции в мотиве допускаются несколько букв. Например, A/G - две буквы, A/G/T - три буквы. Для задания мотивов используется входной файл. Для кодировки комбинаций аминокислот используется 15 буквенный IUPAC код (см. статью MainText.pdf в материалах).  Пример файла мотивов находится в материалах - "FoxA2_major_30.mot".

Необходимо посчитать количество совпадений мотива и частоту (кол-во совпадений/ кол-во массивов) для данного ChIP-seq (одно или несколько совпадения в одном массиве приравнивается к совпадению, т.е. достаточно найти первое совпадения мотива с частью массива).

## Требования

- C++23
- CMake 3.20+
- MPI (OpenMPI, MPICH, или Intel MPI)
- OpenMP
- Google Test (для тестирования)

## Сборка

```bash
# Клонирование репозитория
git clone <https://github.com/vladimir-ponomarenko/SibsutisParallelProgramOptimization.git>
cd SibsutisParallelProgramOptimization

# Сборка проекта
./build.sh
```

## Использование

```bash
# Запуск с MPI
mpirun -n 4 ./DNAMotifFinder sequences.fst motifs.mot

# Сохранение результатов в файл
mpirun -n 4 ./DNAMotifFinder sequences.fst motifs.mot results.txt

# Настройка количества OpenMP потоков
mpirun -n 4 ./DNAMotifFinder --threads 8 sequences.fst motifs.mot

# Подробный вывод
mpirun -n 4 ./DNAMotifFinder --verbose sequences.fst motifs.mot
```

### Параметры командной строки

- `-h, --help` - Показать справку
- `-t, --threads <num>` - Количество OpenMP потоков на процесс
- `-v, --verbose` - Вывод с статистикой производительности

### Формат входных файлов

#### ChIP-seq файл (.fst)
```
>sequence_id	metadata1	metadata2	...
SEQUENCE_LINE_1
SEQUENCE_LINE_2
...
```

#### Файл мотивов (.mot)
```
MOTIF_PATTERN	score1	score2	score3
TRTWKACH	66.9	15.278	1462.3
RTTKACHY	54.3	19.8811	620.023
...
```

## Архитектура

### Основные компоненты

1. **IUPACCodes** - Управление IUPAC кодами для неоднозначных нуклеотидов
2. **DNAParser** - Парсинг входных файлов
3. **MotifFinder** - Основной алгоритм поиска мотивов
4. **MPIManager** - Управление MPI коммуникацией
5. **ParallelProcessor** - Координация MPI и OpenMP


## Производительность

### Параллелизация

- **MPI**: Распределение последовательностей между процессами
- **OpenMP**: Параллельная обработка мотивов внутри процесса
- **Оптимизация**: Минимизация коммуникации и эффективное использование кэша

## Тестирование

```bash
# Запуск всех тестов
make test

# Запуск конкретных тестов
./dna_motif_tests --gtest_filter=IUPACCodesTest.*
./dna_motif_tests --gtest_filter=DNAParserTest.*
./dna_motif_tests --gtest_filter=MotifFinderTest.*

# MPI тесты
mpirun -n 2 ./dna_motif_tests --gtest_filter=MPIManagerTest.*
mpirun -n 2 ./dna_motif_tests --gtest_filter=ParallelProcessorTest.*
```

## Примеры

### Пример 1: Базовый поиск мотивов

```bash
mpirun -n 4 ./DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot
```

### Пример 2: Оптимизированная конфигурация

```bash
mpirun -n 8 ./DNAMotifFinder --threads 4 --verbose data/FoxA2_5000.fst data/FoxA2_major_30.mot results.txt
```

## Результаты

Приложение выводит таблицу с результатами:
- **Motif Pattern** - Паттерн мотива
- **Match Count** - Количество последовательностей с совпадениями
- **Frequency** - Частота совпадений (0.0 - 1.0)

```bash
Motif_Pattern	Match_Count	Frequency
TRTWKACH    	942	        0.376800
RTTKACHY	    738	        0.295200
TMAAYANS	    730	        0.292000
TWKACHYW	    781	        0.312400
TTKRTYTW	    417	        0.166800
TKAHYTWK	    642	        0.256800
TRTTKRTY	    377	        0.150800
RAAYHAAY	    428	        0.171200
```