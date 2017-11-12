function lab_1()
    
    DEBUG = 1;
    
    MAXIMIZE = 1;
    
    C = dlmread('data.txt');
    
    [m, n] = size(C);
    
    if (m ~= n)
        error('The matrix must be square')
    end

    n = size(C, 1);    
    Cm = C;
    
    % Координаты СНН
    independent_zeros = zeros(1, n);        
    
    fprintf('Исходная матрица стоимостей:\n');
    print_Cm(C, independent_zeros);
    
    % В случае задачи максимизации вычитаем элементы из максимума столбца
    if MAXIMIZE ~= 0        
        for i = 1:n
            column_max = max(Cm(i, :));
            for j = 1:n
                Cm(i, j) = column_max - Cm(i, j);
            end
        end
        if DEBUG ~= 0
            fprintf('Переход к задаче максимизации\n');
            print_Cm(Cm, independent_zeros);
        end
    end
            
    % В каждой строке вычитаем минимум этой строки
    for i = 1:n
        minimal = min(Cm(i, :));
        for j = 1:n
            Cm(i, j) = Cm(i, j) - minimal;
        end
    end
    if DEBUG ~= 0
        fprintf('Вычитание минимума строки в каждой строке\n');
        print_Cm(Cm, independent_zeros);
    end    
    
    % В каждом столбце вычитаем миниум этого столбца
    for i = 1:n
        column_min = min(Cm(:, i));
        for j = 1:n
            Cm(j, i) = Cm(j, i) - column_min;
        end
    end
    if DEBUG ~= 0
        fprintf('Вычитание минимума столбца в каждом столбце\n');
        print_Cm(Cm, independent_zeros);
    end     
    
    % Строим начальную СНН
    for i = 1:n
        j = 1;
        flag = 0;
        while j <= n && flag == 0
            if Cm(j, i) == 0 && sum(independent_zeros == j) == 0
                independent_zeros(i) = j;
                flag = 1;
            end
            j = j + 1;
        end
    end
    if DEBUG ~= 0
        fprintf('Выбор независимых нулей\n');
        print_Cm(Cm, independent_zeros);
    end      
   
    % Выделенные столбцы
    marked_columns = zeros(1, n);
    
    % Выделяем столбцы, содержащие независимые нули
    for i = 1:n
        if independent_zeros(i) > 0
            marked_columns(i) = 1;
        end
    end
    
    % Выделенные строки
    marked_strings = zeros(1, n);
        
    % Координаты штрихованных нулей
    marked_zeros = zeros(1, n);    
    
    iter = 1;
    if DEBUG ~= 0
        fprintf('Итерация %d: отметка столбцов\n', iter);
        print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
        fprintf('Колличество помеченных столбцов k = %d\n', sum(independent_zeros > 0));
    end 
    
    % Цикл до тех пор, пока решение не будет оптимальным
    while sum(independent_zeros == 0) ~= 0
        
        % Строка обнаруженного невыделенного нуля
        new_zero_string = 0;
        
        % Отмечаем невыделенный нуль штрихом
        i = 1;
        flag = 0;
        for i = 1 : n
            if marked_columns(i) == 0
                for j = 1 : n
                    if Cm(j, i) == 0 && marked_strings(j) == 0
                        new_zero_string = j;
                        marked_zeros(i) = j;
                        flag = 1;
                        break;
                    end
                end
            end
            if flag == 1
                break;
            end
        end
        
        % Обнаружили невыделенный нуль
        if flag == 1
            
            % Проверяем, что в строке обнаруженного нуля нет 0*
            if sum(independent_zeros == new_zero_string) == 0
                
                % Строим L-цепочку
                l_columns = find(marked_zeros > 0);
                l_chain = zeros(size(l_columns, 2) * 2 - 1, 2);
                l_chain(1, :) = [marked_zeros(l_columns(1)) l_columns(1)];
                for i = 2:size(l_columns, 2)
                    l_chain(i * 2 - 2, :) = [independent_zeros(l_columns(i - 1)) l_columns(i - 1)];
                    l_chain(i * 2 - 1, :) = [marked_zeros(l_columns(i)) l_columns(i)];
                end
                        
                % Удаляем штрихи и отметки строк и столбцов
                marked_zeros = zeros(1, n);
                marked_strings = zeros(1, n);
                marked_columns = zeros(1, n);                
                
                % Формируем новую СНН
                for i = 1:size(l_chain, 1)
                    if mod(i, 2) ~= 0
                        independent_zeros(l_chain(i, 2)) = l_chain(i, 1);
                    end
                end
                
                % Заново выделяем столбцы
                for i = 1:n
                    if independent_zeros(i) > 0
                        marked_columns(i) = 1;
                    end
                end     
                
                if DEBUG ~= 0                  
                    fprintf('-----------------------\n');
                    fprintf('Итерация %d: новые независимые нули\n', iter);
                    print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
                    fprintf('Колличество помеченных столбцов k = %d\n', sum(independent_zeros > 0));                    
                end
                iter = iter + 1;
            else
                % Переносим отметку со столбца на строку
                marked_strings(new_zero_string) = 1;
                marked_columns(independent_zeros == new_zero_string) = 0;
                if DEBUG ~= 0                  
                    fprintf('Итерация %d: перенос отметки\n', iter);
                    print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
                end
            end
        else
            % Невыделенных нулей необнаружено, создаем новые нули
            
            minimal = Inf();
            for i = 1:n
                for j = 1:n
                    if marked_strings(i) == 0 && marked_columns(j) == 0 && Cm(i, j) < minimal
                        minimal = Cm(i, j);
                    end
                end
            end
            for i = 1:n
                for j = 1:n
                    if marked_strings(i) == 1
                        Cm(i, j) = Cm(i, j) + minimal;
                    end
                    if marked_columns(j) == 0
                        Cm(i, j) = Cm(i, j) - minimal;
                    end
                end            
            end
            
            if DEBUG ~= 0
                fprintf('Минимальный элемент h = %d\n', minimal);
                fprintf('Итерация %d: после изменения элементов\n', iter);
                print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
            end
            
            % Обновляем СНН
%             for i = 1:n
%                 j = 1;
%                 flag = 0;
%                 while j <= n && flag == 0
%                     if Cm(j, i) == 0 && sum(independent_zeros == j) == 0
%                         independent_zeros(i) = j;
%                         flag = 1;
%                     end
%                     j = j + 1;
%                 end
%             end
        end
    end
    
    fprintf('-----------------------\n');
    fprintf('Матрица назначений\n');
    % Печать матрицы назначений
    for i = 1:n
        for j = 1:n
            if independent_zeros(j) == i
                fprintf('    1  ');
            else
                fprintf('    0  ');
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
    
    % Вывод итогового значения
    s = zeros(1, n);
    for i = 1:n
        if (independent_zeros(i) > 0)
            s(i) = C(independent_zeros(i), i);
        end
    end
    
    fprintf('\nНазначенные работы: ');
    disp(s);
    
    if MAXIMIZE ~= 0
        fprintf('    Max: %d\n\n', sum(s));
    else
        fprintf('   Min: %d\n\n', sum(s));
    end
end
 
% Печать матрицы стоимостей
function print_Cm(Cm, independent_zeros)
    n = size(Cm, 1);
    for i = 1:n
        for j = 1:n
            if independent_zeros(j) == i
                fprintf('  %3d* ', Cm(i, j));
            else
                fprintf('  %3d  ', Cm(i, j));
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end
    
function [] = print_matrix_full(mat, asn, mcols, mrows, mzeros)
    [~, n] = size(mat);
    for col = 1 : n
        if (mcols(col) == 1)
            fprintf('  +  ');
        else
            fprintf('     ');
        end
    end
    fprintf('\n'); 
    for row = 1 : n
        for col = 1 : n
            fprintf('%3d', mat(row, col));
            if (asn(col) == row)
                fprintf('* ');
            else
                if (mzeros(col) == row)
                    fprintf('` ');
                else
                    fprintf('  ');
                end
            end
        end

        if (mrows(row) == 1)
            fprintf(' +');
        end

        fprintf('\n');
    end
end







