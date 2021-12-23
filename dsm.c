#include <mpi.h>
#include <stdio.h>
//Применение модификации к переменным
void operation(int* mod, int* vars){	
	switch(mod[1]){
	case 0:
		vars[mod[0]]+=mod[2];
		break;
	case 1:
		vars[mod[0]]-=mod[2];
		break;
	case 2:
		vars[mod[0]]*=mod[2];
		break;
	case 3:
		vars[mod[0]]/=mod[2];
		break;
	case 4:
		vars[mod[0]]%=mod[2];
		break;
	default:
		printf("cant use it");
		break;
	}
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    
    //Инициализация модифицируемых переменных
	int vars[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	//Структура модификации
	int mod[3] = {0,0,0};
	
	//Структура модификации автора
	int cur_mod[3] = {0,0,0};
	
	//Количество процессов
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    //Служебные переменные
    int changes = 1; //сколько изменений делает каждый процесс
    int rank;
    int change_num = 0;
    int cur_num = 0;
    int ask = 1;
    MPI_Status status;
    MPI_Request request;
    
    //Получаем ранк процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //Процесс координатор
	if (rank == 0){
		for(int i = 0; i < (world_size-1) * changes; i++){
			//Принимаем запрос на модификацию
			MPI_Recv (&ask, 1, MPI_INT, MPI_ANY_SOURCE, 403, MPI_COMM_WORLD, &status);
			//Присваиваем номер запросу и отправляем автору	
			change_num += 1;
			MPI_Send(&change_num, 1, MPI_INT, status.MPI_SOURCE, 401, MPI_COMM_WORLD);
			//Получаем от автора модификацию
			MPI_Recv (mod, 3, MPI_INT, MPI_ANY_SOURCE, 402, MPI_COMM_WORLD, &status);
			//Рассылаем модификацию всем (в том числе и автору)
			for (int j = 1; j < world_size; j++){
				MPI_Send(mod, 3, MPI_INT, j, change_num, MPI_COMM_WORLD);		
			}
			//Принимаем актуальную версию переменных
			MPI_Recv (vars, 10, MPI_INT, status.MPI_SOURCE, 404, MPI_COMM_WORLD, &status);
		}
		//Модификация на координирующем процессе
		change_num += 1;
		cur_mod[1] = 0;
		cur_mod[0] = rank;
		cur_mod[2] = 2;
		operation(cur_mod, vars);
		for (int j = 1; j < world_size; j++){
				MPI_Send(cur_mod, 3, MPI_INT, j, change_num, MPI_COMM_WORLD);		
		}
		//Рассылаем всем процессам, что модификации закончились
		for (int j = 1; j < world_size; j++){
			MPI_Send(&change_num, 1, MPI_INT, j, 777, MPI_COMM_WORLD);		
		}
		//Выводим последнюю версию переменных.
	    for(int i = 0; i < 10; i++){
			printf("%d ", vars[i]);
		}
		printf(" - final version.\n");
		//Получаем от всех процессов их версии переменных (должны совпадать)
	    for(int i = 1; i < world_size; i++){
			MPI_Recv (vars, 10, MPI_INT, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status);
			for(int i = 0; i < 10; i++){
				printf("%d ", vars[i]);
			}
			printf(" - version from %d process!\n", status.MPI_SOURCE);
		}
		
	//Остальные процессы
	} else {
		for(int i = 0; i < changes; i++){
			//Какую модификацию применяем
			//0 - сумма
			//1 - разность
			//2 - умножение
			//3 - целочисленное деление
			//4 - взятие остатка
			cur_mod[1] = 0;
			//К какой переменной
			cur_mod[0] = rank;
			//Число на которое изменяется (прибавить это число, умножить на это число и тд)
			cur_mod[2] = 2;
			// Отправляем запрос на модификацию
			MPI_Send(&ask, 1, MPI_INT, 0, 403, MPI_COMM_WORLD);
			// Принимаем порядковый номер нашей модификации
			MPI_Recv (&change_num, 1, MPI_INT, 0, 401, MPI_COMM_WORLD, &status);
			// Принимаем все модификации до нашей
			for (int j = cur_num+1; j < change_num; j++){
				MPI_Recv (mod, 3, MPI_INT, 0, j, MPI_COMM_WORLD, &status);
				operation(mod, vars);
			}
			// Теперь мы считаем, что у нас версия за один до нашей модификации
			cur_num = change_num - 1;
			// Отсылаем нашу модификацию
			MPI_Send(cur_mod, 3, MPI_INT, 0, 402, MPI_COMM_WORLD);
			// Принимаем любую модификацию (такая будет точно как минимум одна, так как мы только что отправили свою)
			MPI_Recv (cur_mod, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			// Принимаем все модификации до пришедшей (в теге указан ее номер, может быть как наша, так и уже более поздняя)
			for (int j = cur_num+1; j < status.MPI_TAG; j++){
				MPI_Recv (mod, 3, MPI_INT, 0, j, MPI_COMM_WORLD, &status);
				operation(mod, vars);
			}
			// Применяем пришедшую модификацию
			operation(cur_mod, vars);
			// Теперь мы считаем, что у нас последняя актуальная версия (>= нашей модификации)
			cur_num = status.MPI_TAG;
			//Отправляем последнюю версию координатору
			MPI_Send(vars, 10, MPI_INT, 0, 404, MPI_COMM_WORLD);
			
		}
		//Получаем от координатора извещение, что модификации закончились и номер какая сейчас самая актуальная версия
		MPI_Recv (&change_num, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, &status);
		//Если в процессе не самая актуальная, обновляемся до нее
		for (int j = cur_num+1; j <= change_num; j++){
			MPI_Recv (mod, 3, MPI_INT, 0, j, MPI_COMM_WORLD, &status);
			operation(mod, vars);
		}
		//Отправляем координатору итоговую версию параметров (у всех процессов должен быть одинаковый итог)
		MPI_Send(vars, 10, MPI_INT, 0, 4, MPI_COMM_WORLD);
	}
	
    MPI_Finalize();
}
