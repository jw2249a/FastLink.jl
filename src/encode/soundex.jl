module Soundex
export soundex
import StaticStrings: StaticString

function encode(chr::Char)::Char
    if in(chr, ['a', 'e', 'i', 'o','y', 'u'])
        return '9'
    elseif in(chr,['s','c','k','j','g','z','x','q'])
        return '2'
    elseif in(chr,['n','m'])
        return '5'
    elseif 'l' == chr
        return '4'
    elseif 'r' == chr
        return '6'
    elseif in(chr,['t','d'])
        return '3'
    elseif in(chr,['h','w'])
        return '8'
    elseif in(chr,['v','b','p','f'])
        return '1'
    elseif chr === ' '
        return '7'
    else
        @error "Unknown character encountered $chr"
    end
end

function soundex(s::T)::StaticString{4} where T <: AbstractString
    output = ['0' for _ in 1:4]
    index = 1
    chr=lowercase(s[1])
    previous_encoding = encode(chr)
    output[index] = uppercase(chr)
    for chr in s[2:end]
        if index <= 3
            encoding = encode(lowercase(chr))
            if encoding === '7'
                # continue if space
                continue
            end
            # if vowel or ignorable vars 'h' or 'w'
            if encoding !== '9' && encoding !== '8'
                if encoding !== previous_encoding
                    # Rule 4 on Consonant Separators
                    if encoding === output[index]
                        # If a vowel (A, E, I, O, U) separates two consonants that have the same code, code consonant to the right of vowel
                        if previous_encoding === '9'
                            index += 1
                            output[index] = encoding
                        end
                        # If "H" or "W" separate two consonants that have code, the consonant to the right of the vowel is not coded
                    else
                        # different letter code soundex
                        index += 1
                        output[index] = encoding
                    end
                end
            end
        else
            break
        end
        previous_encoding = encoding
    end
    return StaticString(String(output))
end
end # module
