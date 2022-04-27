%filename HH.m
%numerical solution of the space-clamped Hodgkin-Huxley equations
clear all
clf
global check;
global t1p t2p ip; %parameters for the function izero(t)
in_HH
in_mhnv

%INPUTS
input1 = 1; %and, or, not, xor
input2 = 0; %and, or, xor
circuit = 3; %1=and, 2=or, 3=not, 4=xor

andNum = 4;
orNum = 4;
notNum = 3;

andArr = [0 0 0; 0 0 0; 0 0 0; 1 2 3];
orArr = [0 0 0; 0 0 0; 0 0 0; 1 2 3];
notArr = [0 0 0; 0 0 0; 1 2 0];

andIp = [0, 9*input1, 9*input2];
orIp = [10, 9*input1, 9*input2];
notIp = [20, 9*input1];

boolNum = 12;
boolArr = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0;
           3 4 0; 3 5 0; 8 7 1; 9 6 1; 10 11 2];
boolIp = [0, 10, 20, 9*input1, 9*input2, 9*input1, 9*input2];

inhib = 10;
inhib2 = 10;
wait = -1;
wait2 = -1;

if (circuit == 1)
    numNeurons = andNum;
    ipArr = andIp;
    connectArr = andArr;
elseif (circuit == 2)
    numNeurons = orNum;
    ipArr = orIp;
    connectArr = orArr;
elseif (circuit == 3)
    numNeurons = notNum;
    ipArr = notIp;
    connectArr = notArr;
    inhib = 2;
else
    numNeurons = boolNum;
    ipArr = boolIp;
    connectArr = boolArr;
    inhib = 4;
    inhib2 = 5;
    wait = 6;
    wait2 = 2;
end



cArr = ones(numNeurons) * c;
vArr = ones(numNeurons) * vstart; 
gSyn = 0;

mArr = ones(numNeurons) * m;
hArr = ones(numNeurons) * h;
nArr = ones(numNeurons) * n;

m_plot=zeros(numNeurons,klokmax);
h_plot=zeros(numNeurons,klokmax);
n_plot=zeros(numNeurons,klokmax);
v_plot=zeros(numNeurons,klokmax);
t_plot=zeros(1,klokmax);
c_plot=zeros(numNeurons,klokmax);


for klok=1:klokmax
    
    for curr=1:numNeurons
      t=klok*dt;                      %note time
      mArr(curr)=snew(mArr(curr),alpham(vArr(curr)),betam(vArr(curr)),dt); %update m
      hArr(curr)=snew(hArr(curr),alphah(vArr(curr)),betah(vArr(curr)),dt); %update h
      nArr(curr)=snew(nArr(curr),alphan(vArr(curr)),betan(vArr(curr)),dt); %update n
      gNa=gNabar*(mArr(curr)^3)*hArr(curr);    %sodium conductance
      gK =gKbar*(nArr(curr)^4);    %potassium conductance
      
      if connectArr(curr, 1) == 0
        g=gNa+gK+gLbar;         %total conductance
        gE=gNa*ENa+gK*EK+gLbar*EL;       %gE=g*E
        cArr(curr) = snew(cArr(curr), alphac(vArr(curr)), betac(vArr(curr)), dt);
      else
        gSum = 0;
        eSum = 0;
        for cc=1:3
            if (connectArr(curr,cc) ~= 0)
                temp = gSynbar * cArr(connectArr(curr, cc))^5 * ((2*vArr(connectArr(curr, cc))/25)/(exp(2*vArr(connectArr(curr, cc))/25)-1));
                gSum = gSum + temp;
                if connectArr(curr,cc) == inhib || connectArr(curr,cc) == inhib2
                    eSum = eSum + (temp * -100);
                else
                    eSum = eSum + (temp * 0);
                end
            end
        end
        g = gNa + gK + gLbar + gSum;
        gE=gNa*ENa+gK*EK+gLbar*EL+ eSum;
        cArr(curr) = snew(cArr(curr), alphac(vArr(curr)), betac(vArr(curr)), dt);
      end
      
      %post synaptic, one for each going in (gSyn)
      
      %save old value of v for checking purposes:
      v_old=vArr(curr);
      %update v:
      if (connectArr(curr, 1) == 0)
        ip = ipArr(curr)
      end
      
      if (curr == wait || curr == wait+1) 
          if (t > 6)
              vArr(curr)=(vArr(curr)+(dt/C)*(gE+ip))/(1+(dt/C)*g);
          end
      elseif (curr == wait2) 
          if (t > 8)
              vArr(curr)=(vArr(curr)+(dt/C)*(gE+ip))/(1+(dt/C)*g);
          end
      elseif curr == 1 || connectArr(curr, 1) == 0
          vArr(curr)=(vArr(curr)+(dt/C)*(gE+ip))/(1+(dt/C)*g);
      else
          vArr(curr)=(vArr(curr)+(dt/C)*(gE))/(1+(dt/C)*g);  
      end

      if(check)
        E=gE/g;
        chv=C*(vArr(curr)-v_old)/dt+g*(v-E)-izero(t);
      end
      %store results for future plotting: 
           
      if curr==1
          t_plot(klok)=t;
      end
      v_plot(curr,klok)=vArr(curr);
      c_plot(curr,klok)=cArr(curr);
      m_plot(curr,klok)=mArr(curr);
      h_plot(curr,klok)=hArr(curr);
      n_plot(curr,klok)=nArr(curr);
      
    end
end


if numNeurons > 10
    for cn= 2:numNeurons 
        titles = ["AK", "OK", "NK", "A", "B", "A'", "B'", "NA", "NB", "T1", "T2", "O"];
        subplot(11,1,cn-1),plot(t_plot,v_plot(cn,:),'color',[0,(cn/numNeurons)*1,1], 'LineWidth',1)
        title(titles(cn));

        ylim([-100, 50]);
        hold on
    end
else
    for cn= 1:numNeurons  
        subplot(4,1,cn),plot(t_plot,v_plot(cn,:),'color',[0,(cn/numNeurons)*1,1], 'LineWidth',1)
        if cn==1
            title('Clock Neuron');
        elseif cn==numNeurons
            title('Output Neuron');
        else
            title(['Input Neuron ', num2str(cn-1)])
        end
        ylim([-100, 50]);
        hold on
%     subplot(3,1,2),plot(t_plot,c_plot(cn,:),'color',[0,(cn/numNeurons)*1,1])
%     xlabel('Calcium') 
%     subplot(3,1,3),plot(t_plot,m_plot(cn,:),t_plot,h_plot(cn,:),t_plot,n_plot(cn,:)) 
%     xlabel('M,H,N') 
    end
end
