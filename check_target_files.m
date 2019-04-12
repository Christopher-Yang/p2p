file = dlmread(['jump_tgt_files/subj_11/B8.tgt'])
num_jumps = 0

will_jump = file(:,8)'
for n = 1:length(will_jump)
    if will_jump(1,n) == 3
        num_jumps = num_jumps+1
    end
end

right_jumps = 0
left_jumps = 0
jump_direction = file(:,5)'
for n = 1:length(jump_direction)
    if num_jumps ~= 0
        if jump_direction(1,n) == 0.03
            right_jumps = right_jumps + 1
        else
            if jump_direction(1,n) == -0.03
                left_jumps = left_jumps + 1
            end
        end
    end
end

target_x = file(:,2)'
target_y = file(:,3)'
for n = 1:12
    figure(n)
    scatter(target_x(1,(((n-1)*8)+1):(n*8)),target_y(1,(((n-1)*8)+1):(n*8)))
end