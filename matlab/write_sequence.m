function write_sequence(filename, seq, nO, nU)
% write the given io sequence (a 1 x size matrix) to file in tom format
len = size(seq,2);
if nU != 0
    len = len / 2;
end
fid = fopen(filename, 'w');
fprintf(fid, '{"Type":"SEQUENCE",');
fprintf(fid, '"length":%i,', len);
fprintf(fid, '"nU":%i,', nU);
fprintf(fid, '"nO":%i,"data":[', nO);
for i = 1:size(seq, 2)-1
    fprintf(fid, '%i,', seq(1,i));
end
fprintf(fid, '%i]}', seq(1,size(seq,2)));
fclose(fid);

