local str = "hello\
world"
local x = true and false or true
--print(x) -- true

local y = 1 > 2 and true or false
--print(y) -- true
local t1 = { 3, 4, 7, 1, 5, 2, 8, 12, 15, 5}
function selected(t, k)
    local i, temp = 0, 0
    math.randomseed(os.time())
    for j = 1, k do
        i = math.random(1, #t)
        temp = t[i]
        t[i] = t[j]
        t[j] = temp
    end
    local result = {}
    for i = 1, k do
        table.insert(result, t[i])
    end
    return result
end
local t2 = selected(t1, 10)
--[[
for i = 1, #t2 do
    print(t2[i])
end
]]

function insertSort(t, n)
    for i=2, n do
        local key = t[i]
        local j = i - 1
        while j >= 1 and t[j] > key do
            t[j+1] = t[j]
            j = j - 1
        end
        t[j+1]=key
    end
end

function bubbleSort(t, n)
    for i=1, n do
        for j=n, i+1, -1 do
            if t[j-1]>t[j] then
                local temp = t[j]
                t[j] = t[j-1]
                t[j-1] = temp
            end
        end
    end
end

function selecteSort(t, n)
    for i=1, n-1 do
        local k = i
        for j=i+1, n do
            if t[j] < t[k] then
                k = j
            end
        end
        if k ~= i then
            local temp = t[i]
            t[i] = t[k]
            t[k] = temp
        end
    end
end

function partition(t, l, r)
    local k = t[l]
    while l < r  do
        while l < r and t[r] >= k do
            r = r - 1
        end
        t[l] = t[r]
        while l < r and t[l] <= k do
            l = l + 1
        end
        t[r] = t[l]
    end
    t[l] = k
    return l
end

function quickSort(t, l, r)
    if l >= r then return end
    local p = partition(t, l, r)
    quickSort(t, l, p-1)
    quickSort(t, p+1, r)
end

function merge(t, l, m, r)
    if l >= r then return end
    local l1 = m - l + 1
    local l2 = r - m
    local t1, t2 = {}, {}
    for n=1, l1 do
        t1[n] = t[l+n-1]
    end
    for n=1, l2 do
        t2[n] = t[m+n]
    end
    local i, j, k = 1, 1, l
    while i <= l1 and j <= l2 do
        if t1[i] <= t2[j] then 
            t[k] = t1[i]
            i = i + 1
        else
            t[k] = t2[j]
            j = j + 1
        end
        k = k + 1
    end
    if i > l1 then
        while j <= l2 do
            t[k] = t2[j]
            k = k + 1
            j = j + 1
        end
    end

    if j > l2 then
        while i <= l1 do
            t[k] = t1[i]
            k = k + 1
            i = i + 1
        end
    end
end

function mergeSort(t, l, r)
    if l >= r then return end
    local m = math.floor((l + r) / 2)
    mergeSort(t, l, m)
    mergeSort(t, m+1, r)
    merge(t, l, m, r)
end

function heapSort(t, n)
    local heap = {}
    function heap:new(t, n)
        heap.length = n
        heap.size = n
        heap.data = {}
        for i = 1, n do
            self.data[i] = t[i]
        end
    end

    function heap:parent(i)
        return math.floor(i / 2)
    end

    function heap:left(i)
        return 2 * i
    end

    function heap:right(i)
        return 2 * i + 1
    end

    function heap:maxHeap(i)
        local l = self:left(i)
        local r = self:right(i)
        local largest = i
        if l <= self.size and self.data[l] > self.data[largest] then
            largest = l
        end
        if r <= self.size and self.data[r] > self.data[largest] then
            largest = r
        end
        if i ~= largest then
            local temp = self.data[i]
            self.data[i] = self.data[largest]
            self.data[largest] = temp
            self:maxHeap(largest)
        end
    end

    function heap:buildHeap()
        for i = math.floor(self.length / 2), 1, -1 do
            self:maxHeap(i)
        end
    end

    function heap:heapSort()
        for i = self.length, 1, -1 do
            local temp = self.data[i]
            self.data[i] = self.data[1]
            self.data[1] = temp
            self.size = self.size - 1
            self:maxHeap(1)
        end
    end

    heap:new(t, n)
    heap:buildHeap()
    heap:heapSort()

    for i = 1, n do
        t[i] = heap.data[i]
    end
end

--insertSort(t2, 10)
--bubbleSort(t2, 10)
--selecteSort(t2, 10)
--quickSort(t2, 1, 10)
--mergeSort(t2, 1, 10)
heapSort(t2, 10)

for i=1, #t2 do
    print(t2[i])
end