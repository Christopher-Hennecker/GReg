function [mult,Gregs] = analyze_seq(seq,GQs,scores,L,startI)
%analyze_seq analyzes potential G4CRs, and stores the information on each
%G4CR including: multiplicity, starting index, ending index, length,
%#tandem G4s


%first we define our variables
startInd = [];
endGQs = [];
sGQ = [];
mult = zeros(size(seq));


%Next we sequentially go through the G4matrix (GQs) by the length using a
%and check if each G4 motif is present
for i = 1:length(GQs)
    
    %If the size of the G4 matrix is less than the size of the sequence we
    %use a sliding window to go through the full sequence
    if size(GQs{i},2) < length(seq)
        z = length(seq)-size(GQs{i},2);
        for j = 1:z+1
            %Here the index returns a value of 1 when the G4 motif is
            %present
            ind = sum(GQs{i}.*seq(j:j-1+size(GQs{i},2)),2) == 12;
            if sum(scores{i}(ind)) > 0
                sGQ = [sGQ;scores{i}(ind)];
                mult(j:j-1+size(GQs{i},2)) = mult(j:j-1+size(GQs{i},2)) + sum(1*(GQs{i}(ind,:)==1),1);
                endGQs = [endGQs;zeros(size(scores{i}(ind)))+j-1+size(GQs{i},2)];
                startInd = [startInd;zeros(size(scores{i}(ind)))+j];
            end
        end
    end
    %If the size of the G4 matrix is equal to the size of the sequence we
    %just check the whole sequence for G4 motifs
    if size(GQs{i},2) == length(seq)
        ind = sum(GQs{i}.*seq,2) == 12;
        if sum(scores{i}(ind)) > 0
            mult = mult+sum(1*(GQs{i}(ind,:)==1),1);
            sGQ = [sGQ;scores{i}(ind)];
            startInd = [startInd;zeros(size(scores{i}(ind)))+1];
            endGQs = [endGQs;zeros(size(scores{i}(ind)))+length(seq)];
        end
    end
end

%next we go through and find our G4CR(s), which are defined as overlapping
%G4s. I do this using the starting and ending index, finding places where
%there is a gap between the end of the current G4 motif and the start of
%the next
Gregs = {};
if ~isempty(startInd)
    [s,ind] = sort(startInd);
    e = endGQs(ind);
    sI = s(1);
    eI = e(1);
    for i = 2:length(s)
        if s(i) > eI
            Gregs = [Gregs;{mult(sI:eI)},{sI},{eI}];
            sI = s(i);
            eI = e(i);
        elseif s(i) <= eI && e(i) > eI
            eI = e(i);
        end
    end
    Gregs = [Gregs;{mult(sI:eI)},{sI},{eI}];
end

%Next I determine how many G4s can fold in tandem
for i = 1:size(Gregs,1)
    ind = startInd>=Gregs{i,2} & endGQs <= Gregs{i,3};
    Gregs{i,4} = sGQ(ind);
    e = endGQs(ind);
    s = startInd(ind);
    [e,ind] = sort(e);
    s = s(ind);
    numsf = 0;
    while sum(ind) > 0
        v = min(e);
        ind = s>v;
        s = s(ind);
        e = e(ind);
        numsf = numsf+1;
    end
    Gregs{i,5} = numsf;
end

%Finally I update the starting and ending index of the G4CR
for i = 1:size(Gregs,1)
    Gregs{i,2} = Gregs{i,2}+startI-1;
    Gregs{i,3} = Gregs{i,3}+startI-1;
end
end
