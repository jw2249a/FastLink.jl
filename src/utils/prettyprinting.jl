function center_in_line(str::String, width::Int = 100;pad_char::Char = '#')
    str_length = length(str)
    if str_length >= width
        return str
    end

    padding_total = width - str_length
    left_padding = div(padding_total, 2)
    right_padding = padding_total - left_padding

    return repeat(pad_char, left_padding) * str * repeat(pad_char, right_padding)
end
function pretty_print_number(num::Integer)
    num_str = string(num)
    reversed_num_str = reverse(num_str)
    parts = [reversed_num_str[i:min(i+2, end)] for i in 1:3:length(reversed_num_str)]
    return reverse(join(parts, '_'))
end
