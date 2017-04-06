
load dic_atmos.temp.txt

if exist('state.20151218.glob.nc','file')
    state=rdmnc('state.20151218.glob.nc','iter');
elseif exist('state.20151219.glob.nc','file')
    state=rdmnc('state.20151219.glob.nc','iter');
end

dic_atmos=nan(length(state.iter),3);

for i=1:length(state.iter)
   dic_atmos(i,:)=dic_atmos_temp(find(dic_atmos_temp(:,1)==state.iter(i),1,'first'),:);
end

save dic_atmos.glob.txt dic_atmos -ascii -double