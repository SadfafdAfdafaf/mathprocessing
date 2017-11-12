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
    
    % ���������� ���
    independent_zeros = zeros(1, n);        
    
    fprintf('�������� ������� ����������:\n');
    print_Cm(C, independent_zeros);
    
    % � ������ ������ ������������ �������� �������� �� ��������� �������
    if MAXIMIZE ~= 0        
        for i = 1:n
            column_max = max(Cm(i, :));
            for j = 1:n
                Cm(i, j) = column_max - Cm(i, j);
            end
        end
        if DEBUG ~= 0
            fprintf('������� � ������ ������������\n');
            print_Cm(Cm, independent_zeros);
        end
    end
            
    % � ������ ������ �������� ������� ���� ������
    for i = 1:n
        minimal = min(Cm(i, :));
        for j = 1:n
            Cm(i, j) = Cm(i, j) - minimal;
        end
    end
    if DEBUG ~= 0
        fprintf('��������� �������� ������ � ������ ������\n');
        print_Cm(Cm, independent_zeros);
    end    
    
    % � ������ ������� �������� ������ ����� �������
    for i = 1:n
        column_min = min(Cm(:, i));
        for j = 1:n
            Cm(j, i) = Cm(j, i) - column_min;
        end
    end
    if DEBUG ~= 0
        fprintf('��������� �������� ������� � ������ �������\n');
        print_Cm(Cm, independent_zeros);
    end     
    
    % ������ ��������� ���
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
        fprintf('����� ����������� �����\n');
        print_Cm(Cm, independent_zeros);
    end      
   
    % ���������� �������
    marked_columns = zeros(1, n);
    
    % �������� �������, ���������� ����������� ����
    for i = 1:n
        if independent_zeros(i) > 0
            marked_columns(i) = 1;
        end
    end
    
    % ���������� ������
    marked_strings = zeros(1, n);
        
    % ���������� ������������ �����
    marked_zeros = zeros(1, n);    
    
    iter = 1;
    if DEBUG ~= 0
        fprintf('�������� %d: ������� ��������\n', iter);
        print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
        fprintf('����������� ���������� �������� k = %d\n', sum(independent_zeros > 0));
    end 
    
    % ���� �� ��� ���, ���� ������� �� ����� �����������
    while sum(independent_zeros == 0) ~= 0
        
        % ������ ������������� ������������� ����
        new_zero_string = 0;
        
        % �������� ������������ ���� �������
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
        
        % ���������� ������������ ����
        if flag == 1
            
            % ���������, ��� � ������ ������������� ���� ��� 0*
            if sum(independent_zeros == new_zero_string) == 0
                
                % ������ L-�������
                l_columns = find(marked_zeros > 0);
                l_chain = zeros(size(l_columns, 2) * 2 - 1, 2);
                l_chain(1, :) = [marked_zeros(l_columns(1)) l_columns(1)];
                for i = 2:size(l_columns, 2)
                    l_chain(i * 2 - 2, :) = [independent_zeros(l_columns(i - 1)) l_columns(i - 1)];
                    l_chain(i * 2 - 1, :) = [marked_zeros(l_columns(i)) l_columns(i)];
                end
                        
                % ������� ������ � ������� ����� � ��������
                marked_zeros = zeros(1, n);
                marked_strings = zeros(1, n);
                marked_columns = zeros(1, n);                
                
                % ��������� ����� ���
                for i = 1:size(l_chain, 1)
                    if mod(i, 2) ~= 0
                        independent_zeros(l_chain(i, 2)) = l_chain(i, 1);
                    end
                end
                
                % ������ �������� �������
                for i = 1:n
                    if independent_zeros(i) > 0
                        marked_columns(i) = 1;
                    end
                end     
                
                if DEBUG ~= 0                  
                    fprintf('-----------------------\n');
                    fprintf('�������� %d: ����� ����������� ����\n', iter);
                    print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
                    fprintf('����������� ���������� �������� k = %d\n', sum(independent_zeros > 0));                    
                end
                iter = iter + 1;
            else
                % ��������� ������� �� ������� �� ������
                marked_strings(new_zero_string) = 1;
                marked_columns(independent_zeros == new_zero_string) = 0;
                if DEBUG ~= 0                  
                    fprintf('�������� %d: ������� �������\n', iter);
                    print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
                end
            end
        else
            % ������������ ����� ������������, ������� ����� ����
            
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
                fprintf('����������� ������� h = %d\n', minimal);
                fprintf('�������� %d: ����� ��������� ���������\n', iter);
                print_matrix_full(Cm, independent_zeros, marked_columns, marked_strings, marked_zeros);
            end
            
            % ��������� ���
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
    fprintf('������� ����������\n');
    % ������ ������� ����������
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
    
    % ����� ��������� ��������
    s = zeros(1, n);
    for i = 1:n
        if (independent_zeros(i) > 0)
            s(i) = C(independent_zeros(i), i);
        end
    end
    
    fprintf('\n����������� ������: ');
    disp(s);
    
    if MAXIMIZE ~= 0
        fprintf('    Max: %d\n\n', sum(s));
    else
        fprintf('   Min: %d\n\n', sum(s));
    end
end
 
% ������ ������� ����������
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







