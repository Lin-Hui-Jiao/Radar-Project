// json.hpp — minimal JSON reader/writer (standard library only). Copy of pshc/include/pshc/json.hpp.
#pragma once
#include <cctype>
#include <cmath>
#include <locale>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pshc {
namespace json {

struct Value {
    enum Type { Null, Bool, Number, String, Array, Object } type = Null;
    bool b = false;
    double num = 0;
    std::string str;
    std::vector<Value> arr;
    std::vector<std::pair<std::string, Value>> obj;   // keeps key order

    bool is_null() const { return type == Null; }
    bool has(const std::string& k) const { return find(k) != nullptr; }
    const Value* find(const std::string& k) const {
        if (type != Object) return nullptr;
        for (auto& kv : obj) if (kv.first == k) return &kv.second;
        return nullptr;
    }
    const Value& operator[](const std::string& k) const {
        const Value* v = find(k);
        if (!v) throw std::runtime_error("json: missing key '" + k + "'");
        return *v;
    }
    const Value& operator[](size_t i) const {
        if (type != Array || i >= arr.size()) throw std::runtime_error("json: index out of range");
        return arr[i];
    }
    size_t size() const { return type == Array ? arr.size() : type == Object ? obj.size() : 0; }
    double as_num() const { if (type != Number) throw std::runtime_error("json: not a number"); return num; }
    int64_t as_int() const { return (int64_t)std::llround(as_num()); }
    const std::string& as_str() const { if (type != String) throw std::runtime_error("json: not a string"); return str; }
    bool as_bool() const { if (type != Bool) throw std::runtime_error("json: not a bool"); return b; }
};

class Parser {
public:
    explicit Parser(const std::string& s) : s_(s) {}
    Value parse() {
        Value v = value(); ws();
        if (p_ != s_.size()) fail("trailing characters");
        return v;
    }

private:
    const std::string& s_;
    size_t p_ = 0;
    [[noreturn]] void fail(const char* m) const {
        throw std::runtime_error(std::string("json: ") + m + " at offset " + std::to_string(p_));
    }
    void ws() { while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' || s_[p_] == '\r')) ++p_; }
    bool lit(const char* w) {
        size_t n = std::char_traits<char>::length(w);
        if (s_.compare(p_, n, w) == 0) { p_ += n; return true; }
        return false;
    }
    Value value() {
        ws();
        if (p_ >= s_.size()) fail("unexpected end");
        Value v;
        char c = s_[p_];
        if (c == '{') {
            v.type = Value::Object; ++p_; ws();
            if (p_ < s_.size() && s_[p_] == '}') { ++p_; return v; }
            for (;;) {
                ws(); if (p_ >= s_.size() || s_[p_] != '"') fail("expected key");
                std::string k = string(); ws();
                if (p_ >= s_.size() || s_[p_] != ':') fail("expected ':'");
                ++p_;
                v.obj.emplace_back(std::move(k), value()); ws();
                if (p_ < s_.size() && s_[p_] == ',') { ++p_; continue; }
                if (p_ < s_.size() && s_[p_] == '}') { ++p_; break; }
                fail("expected ',' or '}'");
            }
        } else if (c == '[') {
            v.type = Value::Array; ++p_; ws();
            if (p_ < s_.size() && s_[p_] == ']') { ++p_; return v; }
            for (;;) {
                v.arr.push_back(value()); ws();
                if (p_ < s_.size() && s_[p_] == ',') { ++p_; continue; }
                if (p_ < s_.size() && s_[p_] == ']') { ++p_; break; }
                fail("expected ',' or ']'");
            }
        } else if (c == '"') {
            v.type = Value::String; v.str = string();
        } else if (lit("true")) { v.type = Value::Bool; v.b = true; }
        else if (lit("false")) { v.type = Value::Bool; v.b = false; }
        else if (lit("null")) { v.type = Value::Null; }
        else {
            size_t q = p_;
            while (q < s_.size() && (std::isdigit((unsigned char)s_[q]) || s_[q] == '-' || s_[q] == '+' || s_[q] == '.' ||
                                     s_[q] == 'e' || s_[q] == 'E')) ++q;
            if (q == p_) fail("unexpected character");
            v.type = Value::Number;
            std::istringstream is(s_.substr(p_, q - p_)); is.imbue(std::locale::classic());
            is >> v.num;
            if (is.fail()) fail("bad number");
            p_ = q;
        }
        return v;
    }
    std::string string() {
        ++p_;   // opening quote
        std::string o;
        while (p_ < s_.size() && s_[p_] != '"') {
            char c = s_[p_++];
            if (c != '\\') { o += c; continue; }
            if (p_ >= s_.size()) fail("bad escape");
            char e = s_[p_++];
            switch (e) {
                case '"': o += '"'; break; case '\\': o += '\\'; break; case '/': o += '/'; break;
                case 'b': o += '\b'; break; case 'f': o += '\f'; break; case 'n': o += '\n'; break;
                case 'r': o += '\r'; break; case 't': o += '\t'; break;
                case 'u': {
                    if (p_ + 4 > s_.size()) fail("bad \\u escape");
                    unsigned cp = (unsigned)std::stoul(s_.substr(p_, 4), nullptr, 16); p_ += 4;
                    if (cp < 0x80) o += (char)cp;
                    else if (cp < 0x800) { o += (char)(0xC0 | (cp >> 6)); o += (char)(0x80 | (cp & 0x3F)); }
                    else { o += (char)(0xE0 | (cp >> 12)); o += (char)(0x80 | ((cp >> 6) & 0x3F)); o += (char)(0x80 | (cp & 0x3F)); }
                    break;
                }
                default: fail("bad escape");
            }
        }
        if (p_ >= s_.size()) fail("unterminated string");
        ++p_;
        return o;
    }
};

inline Value parse(const std::string& text) { return Parser(text).parse(); }

inline Value parse_file(const std::string& path) {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) throw std::runtime_error("cannot open " + path);
    std::string s; char buf[65536]; size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, fp)) > 0) s.append(buf, n);
    std::fclose(fp);
    return parse(s);
}

// Streaming writer for reports: w.begin_object(); w.key("a").num(1); ... w.end_object().
class Writer {
public:
    Writer& begin_object() { sep(); os_ << '{'; stack_.push_back(true); first_ = true; ++depth_; return *this; }
    Writer& end_object() { --depth_; nl(); os_ << '}'; stack_.pop_back(); first_ = false; return *this; }
    Writer& begin_array() { sep(); os_ << '['; stack_.push_back(false); first_ = true; ++depth_; return *this; }
    Writer& end_array() { --depth_; os_ << ']'; stack_.pop_back(); first_ = false; return *this; }
    Writer& key(const std::string& k) { sep(); quote(k); os_ << ": "; after_key_ = true; return *this; }
    Writer& num(double v) {
        sep();
        if (!std::isfinite(v)) { os_ << "null"; return *this; }
        char b[40];
        if (v == std::floor(v) && std::fabs(v) < 9.0e15) std::snprintf(b, sizeof b, "%.0f", v);
        else std::snprintf(b, sizeof b, "%.17g", v);
        os_ << b; return *this;
    }
    Writer& num(double v, int sig) {
        sep();
        if (!std::isfinite(v)) { os_ << "null"; return *this; }
        char b[40]; std::snprintf(b, sizeof b, "%.*g", sig, v); os_ << b; return *this;
    }
    Writer& integer(int64_t v) { sep(); os_ << v; return *this; }
    Writer& str(const std::string& s) { sep(); quote(s); return *this; }
    Writer& boolean(bool b) { sep(); os_ << (b ? "true" : "false"); return *this; }
    Writer& null() { sep(); os_ << "null"; return *this; }
    std::string text() const { return os_.str() + "\n"; }
    bool save(const std::string& path) const {
        FILE* fp = std::fopen(path.c_str(), "w"); if (!fp) return false;
        std::string t = text(); std::fwrite(t.data(), 1, t.size(), fp); std::fclose(fp); return true;
    }

private:
    std::ostringstream os_;
    std::vector<bool> stack_;
    bool first_ = true, after_key_ = false;
    int depth_ = 0;
    void nl() { os_ << '\n'; for (int i = 0; i < depth_; ++i) os_ << "  "; }
    void sep() {
        if (after_key_) { after_key_ = false; return; }
        if (stack_.empty()) return;
        const bool was_first = first_;
        if (!first_) os_ << ',';
        first_ = false;
        if (stack_.back()) nl(); else if (!was_first) os_ << ' ';
    }
    void quote(const std::string& s) {
        os_ << '"';
        for (unsigned char c : s) {
            if (c == '"' || c == '\\') os_ << '\\' << c;
            else if (c == '\n') os_ << "\\n";
            else if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); os_ << b; }
            else os_ << c;
        }
        os_ << '"';
    }
};

}  // namespace json
}  // namespace pshc
