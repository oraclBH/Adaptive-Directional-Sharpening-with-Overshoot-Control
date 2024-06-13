function outIM =  adaptiveDirectionalSharpenWithOvershootControl(inIM, win,  gain,  th1,  th2,  strength)

if  ~exist('win', 'var')
    win = 2;
end

if  ~exist('gain', 'var')
    gain = 16;     %【0， 128】
end

if  ~exist('th1', 'var')
    th1 = 10;
end

if  ~exist('th2', 'var')
    th2 = 20;
end
  
if  ~exist('strength', 'var')
    strength = 3;
end

img = double(inIM);
[m, n, s] = size(img);

%% 添加图像扩展函数，使用对称填充
imgpad = zeros(m+2*win, n+2*win, s);
imgpad(win+1:m+win, win+1:n+win, :) = img;
imgpad(1:win, win+1:n+win, :) = img(win+1:-1:2, :, :); 
imgpad(m+1:m+win, win+1:n+win, :) = img(m-1:-1:m-win, :, :); 
imgpad(:,1:win, :) = imgpad(:, 2*win+1:-1:win+2, :);
imgpad(:, n+1:n+win, :) = imgpad(:, n-1:-1:n-win, :);
 
y= 0.299*imgpad(:,:,1) + 0.587*imgpad(:,:,2) + 0.114 *imgpad(:,:,3);
u= - 0.1687*imgpad(:,:,1) - 0.3313*imgpad(:,:,2) + 0.5 *imgpad(:,:,3) + 128;   
v = 0.5*imgpad(:,:,1) - 0.4187*imgpad(:,:,2) - 0.0813 *imgpad(:,:,3) + 128;    
yout = y;
sharpen_map = zeros(size(y));
  
%% 高通滤波器来识别纹理和平坦区域
kernel = [-1,   0,   -1,   0,   -1;
                 0,  -2,   -1,  -2,   0;
                -1,  -1,  20,  -1,  -1;
                 0,   -2,  -1,  -2,   0;
                -1,   0,   -1,   0,  -1];
            
yhp = imfilter(y, kernel, 'replicate');            

%% 使用5*1滤波器
vector_kernel_horiz = [0,  0,  0,  0,  0;
                                     0,  0,  0,  0,  0;
                                    -1,  0,  2,  0, -1;
                                     0,  0,  0,  0,  0;
                                     0,  0,  0,  0,  0;];

vector_kernel_vertic = [ 0,  0,  -1,  0,  0;
                                        0,  0,  0,  0,  0;
                                        0,  0,  2,  0,  0;
                                        0,  0,  0,  0,  0;
                                        0,  0,  -1,  0,  0;];
 

vector_kernel_positiveDiagonal =[-1,  0,  0,  0,  0;
                                                         0,  0,  0,  0,  0;     
                                                         0,  0,  2,  0,  0; 
                                                         0,  0,  0,  0,  0; 
                                                         0,  0,  0,  0,  -1 ];
            
vector_kernel_negativeDiagonal =[ 0,  0,  0,  0,  -1;
                                                         0,  0,  0,  0,  0;     
                                                         0,  0,  2,  0,  0; 
                                                         0,  0,  0,  0,  0; 
                                                        -1,  0,  0,  0,  0];
                                                    
for i = win+1: m+win
    for j = win+1:n+win
        vv = y(i, j);
        vhp = yhp(i, j);
        block = y(i-win:i+win, j-win:j+win);
        vhorz = sum(sum(vector_kernel_horiz.*block));
        vvert = sum(sum(vector_kernel_vertic.*block));
        vpdiag =  sum(sum(vector_kernel_positiveDiagonal.*block));
        vndiag =  sum(sum(vector_kernel_negativeDiagonal.*block));
        ADSOC = [vhorz, vpdiag, vndiag, vvert];
        [maxsharpen, index] = max(ADSOC);
        maxsharpen_gain = ADSOC(5 - index);
        minv = min(block(:));
        maxv = max(block(:));
        alpha = strength / 10;
        
        if vhp > th1
            sharpen_v = maxsharpen_gain * gain / 16;
            if maxsharpen >th2
                tmp = vv + sharpen_v;
                if tmp > maxv
                    tmp = maxv + alpha * (tmp - maxv);
                    sharpen_v = tmp - vv;
                elseif tmp < minv
                    tmp = minv - alpha * (minv - tmp);
                    sharpen_v = tmp - vv;
                end

            end
            sharpen_map(i, j) = sharpen_v;
        end
    end
end

hw = fspecial('average', [3, 3]);
sharpen_map = imfilter(sharpen_map, hw, 'replicate');
yout = y + sharpen_map;
                                                      
y=yout;
R = y + 1.402* (v-128);
G = y - 0.34414* (u-128) - 0.71414* (v-128);
B = y + 1.772 *(u-128);
tmp_img = cat(3,R,G,B);
outIM= tmp_img(win+1:m+win, win+1:n+win, :);
end