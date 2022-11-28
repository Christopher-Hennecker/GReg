function [GQs,Tm_est] = get_GQs(maxL,maxB,minTm)
%get_GQs sequentially creates all possible G4 motifs containing loops of
%length less than maxL, bulges of length less than maxB, and Tm_est greater
%than minTm



%Define outputs: 
%GQs is the cell array of putative G4 motifs
GQs = {};
%G4s is an overflow cell array for GQs.
G4s = {};
T4s = [];
%num is the number of putative G4 motifs.
num = 0;
%Tm_est is a matrix of the estimated melting temperatures
Tm_est = [];


%Sequential creation of putative G4 motifs. This is done using 7 nested for
%loops. 4 for the bulges in the G-tracts (B0/B1/B2/B3) which go up to a
%maximum length (maxB), and incorporates the bulge either between the first
%and second core guanine, or the second and third.
for B0 = 0:maxB*2
    %Creation of the first G-tract - G0
    if B0 == 0
        G0 = [1,1,1];
    elseif B0 > 0 && B0 <=maxB
        G0 = [1,zeros(1,B0)+13,1,1];
    else
        G0 = [1,1,zeros(1,mod(B0,maxB+1)+1)+13,1];
    end
    for L1 = 1:maxL
        %Creation of the first loop - l1
        l1 = zeros(1,L1);
        for B1 = 0:maxB*2
            %Creation of the second G-tract - G1
            if B1 == 0
                G1 = [1,1,1];
            elseif B1 > 0 && B1 <=maxB
                G1 = [1,zeros(1,B1)+13,1,1];
            else
                G1 = [1,1,zeros(1,mod(B1,maxB+1)+1)+13,1];
            end
            for L2 = 1:maxL
                %cCreation of the second loop - l2
                l2 = zeros(1,L2);
                for B2 = 0:maxB*2
                    %Creation of the third G-tract - G2
                    if B2 == 0
                        G2 = [1,1,1];
                    elseif B2 > 0 && B2 <=maxB
                        G2 = [1,zeros(1,B2)+13,1,1];
                    else
                        G2 = [1,1,zeros(1,mod(B2,maxB+1)+1)+13,1];
                    end
                  
                    
                    for L3 = 1:maxL
                        %Creation of the third loop - L3
                        l3 = zeros(1,L3); 
                        for B3 = 0:maxB*2
                            %Creation of the fourth G-tract - G3;
                            if B3 == 0
                                G3 = [1,1,1];
                            elseif B3 > 0 && B3 <=maxB
                                G3 = [1,zeros(1,B3)+13,1,1];
                            else
                                G3 = [1,1,zeros(1,mod(B3,maxB+1)+1)+13,1];
                            end

                            %Calculating the number of total bulged residues
                            B = sum((G0==13)*1)+sum((G1==13)*1)+sum((G2==13)*1)+sum((G3==13)*1);
                            %Calculating the number of G-tracts with bulges
                            tracts = max(G0==13)+max(G1==13)+max(G2==13)+max(G3==13);
                            %Calculating the logarithm of 
                            L = log10(length(l1)*length(l2)*length(l3));
                            %Creating the final G4 motif
                            currGQ = [G0,l1,G1,l2,G2,l3,G3];
                            %Calculating the estimated melting temperature
                            Tm = 89.9-19.2*L-20*tracts-8.5*(B-tracts);
                            if Tm >= minTm
                                %If the estimated temperature is higher
                                %than the minimum Tm, the we store the G4
                                %motif and the Tm
                                Tm_est = [Tm_est;Tm];
                                GQs = [GQs;currGQ];
                                if num > 10000
                                    %if we have stored more than 10,000 G4
                                    %structures we use G4s as an overflow
                                    %cell array to speed up the algorithm
                                    G4s = [G4s;GQs];
                                    T4s = [T4s;Tm_est];
                                    GQs = {};
                                    num = 0;
                                end
                                %to keep track of the number of G4 motifs
                                num = num+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%Combine any overflow G4s with the current GQs
GQs = [GQs;G4s];
T4s = [T4s;Tm_est];
Tm_est = T4s;

%sort all the G4s by length, and update the Tm_est
[~,I] = sort(cellfun(@length,GQs));
Tm_est = Tm_est(I);
GQs = GQs(I);


%Convert the GQs into a sparse array where G4 motifs of the same size are
%combined together in a matrix.
nGQ = {};
nTm_est = {};
currLen = 0;
for i = 1:length(GQs)
    if length(GQs{i}) == currLen
        nGQ{end} = [nGQ{end};sparse(GQs{i})];
        nTm_est{end} = [nTm_est{end};Tm_est(i)];
    else
        nGQ{end+1} = sparse(GQs{i});
        nTm_est{end+1} = Tm_est(i);
        currLen = length(GQs{i});
    end
end


GQs = nGQ;
Tm_est = nTm_est;
end