#!/bin/env bash
#SBATCH --job-name=heat                       # Название задачи
#SBATCH --error=./errs/heat-%j.err            # Файл для вывода ошибок
#SBATCH --output=./logs/heat-%j.log           # Файл для вывода результатов

#SBATCH --account=proj_1593                   # Идентификатор проекта
#SBATCH --time=0:00:30                        # Максимальное время выполнения

#SBATCH --ntasks=4                            # Количество MPI процессов
#SBATCH --ntasks-per-node=4                   # Количество процессов на один узел
#SBATCH --cpus-per-task=1                     # Количество ядер CPU, выделенных для одного процесса
#SBATCH --ntasks-per-socket=4                 # Количество процессов на один сокет
#SBATCH --constraint="type_a|type_b|type_c"   # Предпочитаемый тип узлов
#SBATCH --nodes=1                             # Количество используемых узлов

module load openmpi                           # Загрузка модуля openmpi
mpicxx heat.cpp -o heat
srun ./heat 4 51 0.1
