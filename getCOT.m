%%－－－－－－阶比跟踪算法－－－－－
function [Tn,xtn] = getCOT(x,t,fs,Dmax,pf,order,wu)
    %xtn为输出：等角度采样的信号序列
    %输入：x为等时间间隔采样信号序列，
    %     t为时间，
    %     fs采样频率，
    %     Dmax为最大阶次，
    %     pf为频率曲线序列,
    %     order:拟合频率曲线的阶次，
    %     wu：舍弃的点数

    t = t - min(t);
     
    dw = pi/Dmax; %重采样角度间隔
    dt = 1/fs ;%采样时间间隔
    
    a = polyfit(t,pf,order); %3阶拟合：ft = a(1)*t.^3 + a(2)*t.^2 + a(3)*t+a(4);
    ft = polyval(a,t);  %得到拟合后的频率曲线序列
    Na = length(a);
    for j = 1:Na
        a(j) = a(j)/(Na-j+1);
    end
    
    lenXtn = fix(2*pi*sum(ft*dt)/dw); %计算重采样后的数据长度
    lenXtn = lenXtn -wu;
    
    Tn = zeros(1,lenXtn); % 计算键相时标 
    for ii = 1 : lenXtn      %   求解方程
        temp = ii/(2*Dmax);
        r = roots([a -temp]);
        for kk = 1: length(r)
            if isreal(r(kk)) 
                if r(kk) > 0 && r(kk) < 10000
                    Tn(ii) = r(kk);
                end
            end
        end
    end
    
    xtn = zeros(1,lenXtn);
    
   for ii = 1: lenXtn
       nn = sum(t<Tn(ii));
       xtn(ii) = x(nn) + (x(nn+1)-x(nn)/(t(nn+1)-t(nn)))*(Tn(ii)-t(nn));
   end   
   
%    for ii = 1: length(x)-1
%        for jj = 1:lenXtn
%            if Tn(jj)>=t(ii) && Tn(jj)<=t(ii+1)
%                 xtn(jj) = x(ii) + (x(ii+1)-x(ii)/(t(ii+1)-t(ii)))*(Tn(jj)-t(ii));
%            end
%        end
%    end

end 