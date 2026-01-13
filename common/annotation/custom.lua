function replaceSpacesWithUnderscores(arr)
    if arr and type(arr[1]) == "string" then
        local modifiedString = string.gsub(arr[1], ' ', '_')
        return modifiedString
    else
        return arr and arr[1] or nil
    end
end

function calculateMean(arr)
    if arr == nil then
        print("Error: Input array is nil")
        return 0
    end

    local sum = 0
    local count = 0

    for i, v in ipairs(arr) do
        local value = tonumber(v)
        if value then
            sum = sum + value
            count = count + 1
        end
    end

    if count > 0 then
        return sum / count
    else
        print("Error: No valid numbers found in the array")
        return 0  -- or handle the case when no valid numbers are found
    end
end

function ratio(vals)
    vals = vals[1] -- get 2 values per element. ref and alt counts.
    if vals[2] == 0 then return "0.0" end
    return string.format("%.9f", vals[2] / (vals[1] + vals[2]))
end

function contains(str, tok)
    return string.find(str, tok) ~= nil
end

function split(str, sep)
        local sep, fields = sep or ":", {}
        local pattern = string.format("([^%s]+)", sep)
        str:gsub(pattern, function(c) fields[#fields+1] = c end)
        return fields
end


CLINVAR_SIG = {}
CLINVAR_SIG['0'] = 'uncertain'
CLINVAR_SIG['1'] = 'not-provided'
CLINVAR_SIG['2'] = 'benign'
CLINVAR_SIG['3'] = 'likely-benign'
CLINVAR_SIG['4'] = 'likely-pathogenic'
CLINVAR_SIG['5'] = 'pathogenic'
CLINVAR_SIG['6'] = 'drug-response'
CLINVAR_SIG['7'] = 'histocompatibility'
CLINVAR_SIG['255'] = 'other'
CLINVAR_SIG['.'] = '.'

function intotbl(ud)
    local tbl = {}
    for i=1,#ud do
        tbl[i] = ud[i]
    end
    return tbl
end


function clinvar_sig(vals)
    local t = type(vals)
    -- just a single-value
    if(t == "string" or t == "number") and not contains(vals, "|") then
        return CLINVAR_SIG[vals]
    elseif t ~= "table" then
        if not contains(t, "userdata") then
            if t == "string" then
                vals = split(vals, ",")
            else
                vals = {vals}
            end
        else
            vals = intotbl(vals)
        end
    end
    local ret = {}
    for i=1,#vals do
        if not contains(vals[i], "|") then
            ret[#ret+1] = CLINVAR_SIG[vals[i]]
        else
            local invals = split(vals[i], "|")
            local inret = {}
            for j=1,#invals do inret[#inret+1] = CLINVAR_SIG[invals[j]]
            end
            ret[#ret+1] = table.concat(inret, "|")
        end
    end
    return table.concat(ret, ",")
end
