clear all 
close all 
clc
NB_seq = 20; 
NB_dgt=10; 
NB_states=402; 
NB_sam = 50; 

count=0;
for idx_state=1:2:NB_states-2
   
    if(mod(idx_state,20)==1)
        count=count+1;
        disp(count)
    end
    idx_presented(idx_state) = 25*mod(idx_state, 20)-25+count; 
    idx_presented(idx_state+1) = idx_presented(idx_state);
    
end
idx_presented(NB_states-1:NB_states) = 1; 
plot(idx_presented)

dlmwrite('idx_presented.dat', idx_presented,'-append','delimiter',' ','roffset',1)
