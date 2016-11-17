

count = 0;
repeat = true;
while count < 5 && repeat
    saveOrNot = input('Append pk data to procData? (y/n)','s');
    count = count + 1;
    if strcmpi(saveOrNot,'y')
        procData.Properties.Writable = true;
        procData.hOr = curv_head;
        procData.elicitedSwimInfo = out;
        repeat = false;
    elseif strcmpi(saveOrNot,'n')
        disp('Data not appended to procData!')
        repeat = false;
    else
        disp('Please enter "y" or "n"')
        repeat = true;
    end
end