function str = ternary(condition, trueStr, falseStr)
if condition
    str = trueStr;
else
    str = falseStr;
end
end