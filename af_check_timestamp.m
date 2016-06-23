dt = [];
clear at lt
for n = 1:numel(data);
    dt = [dt, [data(n).Aquadopps(3).Timestamp{:}]-[data(n).Aquadopps(3).LoggerTime{:}]];
    at(n,:) =  [data(n).Aquadopps(3).Timestamp{:}];
    lt(n,:) =  [data(n).Aquadopps(3).LoggerTime{:}];
end
at = at';
lt = lt';
figure(1);clf
subplot(3,1,1)
plot(dt)
subplot(3,1,2)
plot(lt(:))
title('lt')
subplot(3,1,3)
plot(at(:))
title('at')
%%
figure(3);clf
pcolor(at-lt);shading flat
%%
dt = [];
clear at lt
for n = 1:numel(data);
    dt = [dt, [data(n).Microcats(3).Timestamp{:}]-[data(n).Microcats(3).LoggerTime{:}]];
    at(n,:) =  [data(n).Microcats(3).Timestamp{:}];
    lt(n,:) =  [data(n).Microcats(3).LoggerTime{:}];
end
at = at';
lt = lt';
figure(2);clf
subplot(3,1,1)
plot(dt)
subplot(3,1,2)
plot(lt(:))
title('lt')
subplot(3,1,3)
plot(at(:))
title('at')

%%
figure(3);clf
plot(diff(lt(:)))

ddt = diff(lt(:));
ddt=ddt(isfinite(ddt));
unique(ddt)


