function data = check_seq(seq,GQs,scores,L)
seqs = splitSeq(sparse(seq),L);
Gregs = {};
x = zeros(size(seq));
for i = 1:size(seqs,1)
    if sum(seqs{i}) >= 12 
    [x1,gRegsTemp] = analyze_seq(seqs{i,1},GQs,scores,L,seqs{i,2});
    x(seqs{i,2}:seqs{i,3}) = x(seqs{i,2}:seqs{i,3})+x1;
    Gregs = [Gregs;gRegsTemp];
    end
    
end

data{1} = x;
data{2} = Gregs;
end
function seqs = splitSeq(seq,L)
tag = 0;
A = movmean([zeros(1,L+1),seq,zeros(1,L+1)],L+1);

seqs = {};
tag2 = 0;
for j = 1:length(A)

    if A(j) > 0 && tag == 0
        tag = 1;
        startI = j+3;
        tag2 = 1;
    elseif A(j) == 0 && tag == 1
        tag = 0;
        endI = j-5;
        if endI - startI + 1 >= 15 && sum(seq(startI-(L+1):endI-(L+1)))>=12
        seqs = [seqs; {seq(startI-(L+1):endI-(L+1)),startI-(L+1),endI-(L+1)}];
        end
    end
    
end
end


    
