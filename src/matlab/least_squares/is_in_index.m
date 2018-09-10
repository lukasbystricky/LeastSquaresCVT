function result = is_in_index(index_set, v)

result = 0;
for u=index_set.'
    if v == u
        result = 1;
        break
    end
end


end